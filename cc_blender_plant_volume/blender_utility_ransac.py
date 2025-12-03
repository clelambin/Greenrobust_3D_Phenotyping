"""Implementation of the RANSAC algorithm to fit 2D sections for blender process
"""

# Import libraries
from collections.abc import Callable
from dataclasses import dataclass
from statistics import fmean
from copy import deepcopy
from random import choices
from mathutils import Vector, geometry

# Libraries used for test function only
import bpy
import bmesh

# Import user modules
from cc_blender_plant_volume.blender_user_types import Cartesian, Segment

# Script variable used within module
NORMAL_DIR: dict[Cartesian, Cartesian] = {"X":"Y", "Y":"X"}   # Map normal direction of an axis

# Utility functions
def point_to_segment(point:Vector, segment: Segment) -> float:
    """Return the distance from the point to the segment"""
    # Intersect function return the point on the line and the relative normal distance to the line
    intersect, dist = geometry.intersect_point_line(point, *segment)
    # Check if intersected point within the segment by checking if distance is between 0 and 1
    if 0 <= dist <= 1:
        # Intersection point is within the segment
        # Distance is the length of the vector from initial point to the intersect
        return (point-intersect).length
    # Compute direct distance between either point and return minimum
    return min((point-seg).length for seg in segment)

def max_along_axis(points:list[Vector], axis:Cartesian) -> float:
    """Return max value along a given axis from input points"""
    return max(getattr(point, axis.lower()) for point in points)

def min_along_axis(points:list[Vector], axis:Cartesian) -> float:
    """Return min value along a given axis from input points"""
    return min(getattr(point, axis.lower()) for point in points)

def lowest_mid_axis(points:list[Vector], axis:Cartesian) -> float:
    """Return the lowest value along axis for 1st half of points (based on normal axis"""
    normal_axis = NORMAL_DIR[axis]
    # Extract the first half of points along normal axis
    nb_points = len(points)
    points_sorted = sorted(points, key=lambda point:getattr(point, normal_axis.lower()))
    half_points   = points_sorted[0:int(nb_points/2)]
    # Return min value
    return min(getattr(point, axis.lower()) for point in half_points)

def max_normal_axis(points:list[Vector], axis:Cartesian) -> float:
    """Return the axis cooridate for the max along the normal axis"""
    normal_axis = NORMAL_DIR[axis]
    max_point = max(points, key=lambda point:getattr(point, normal_axis.lower()))
    return getattr(max_point, axis.lower())



# Model
@dataclass
class RansacParam:
    """Parameter used to run RANSAC"""
    nb_sample:int = 20
    min_pts_per_line:int = 3
    max_iter:int = 100
    dist_thresh:float = 0.01
    max_fit:float = 0.9


@dataclass
class ModelParam:
    """Parameter with attached value, name and supporting segment"""
    name: str
    direction: Cartesian
    fit_fct: tuple[Callable, Callable]
    values: tuple[float|None, float|None] = (None, None)
    update_fct: Callable = fmean
    bound_min: "ModelParam|None" = None
    bound_max: "ModelParam|None" = None

    def __repr__(self) -> str:
        return f"{self.name} = ({self.values[0]:.3f}, {self.values[1]:.3f})"

    def set_bound(self, limit_min:"ModelParam|None"=None, limit_max:"ModelParam|None"=None) -> None:
        """Set corresponding bounding limits"""
        if limit_min is not None:
            self.bound_min = limit_min
        if limit_max is not None:
            self.bound_max = limit_max

    def supporting_segment(self) -> Segment:
        """Compute segment supporting current parameter"""
        # Extract bounding parameters (to define limits)
        limit_min = 0 if self.bound_min is None else self.bound_min.values[0]
        limit_max = 0 if self.bound_max is None else self.bound_max.values[1]

        # Initialise the points defining the segment
        point_min = Vector((0, 0))
        point_max = Vector((0, 0))

        # Set coordinate along direction to value and normal direction to limit
        normal_axis = NORMAL_DIR[self.direction]
        setattr(point_min, self.direction.lower(), self.values[0])
        setattr(point_max, self.direction.lower(), self.values[1])
        setattr(point_min, normal_axis.lower(), limit_min)
        setattr(point_max, normal_axis.lower(), limit_max)

        # Return segment
        return (point_min, point_max)

    def set_values(self, points:list[Vector]) -> None:
        """Update values by calling each fit_fct function on the input points"""
        new_value = []
        for fct in self.fit_fct:
            new_value.append(fct(points, self.direction))
        self.values = tuple(new_value)


# TODO: Recompute twice for horizontal/vertical parameters (not efficient)
class PotSection:
    """Simplified pot section parameters

    Pot section described by 2 rectangle overlapping each other.
    Each rectangle is described by height (along Y axis) and width (along X axis)
        - Pot : outer rectangle centered in X and offset in Y by ground height
        - Soil : inner rectangle centered in X and with top line overlapping with soil level
    """

    def __init__(self) -> None:
        """Initialise pot section with empty parameter"""

        # Initialise parameters
        pot_width = ModelParam("pot_width", "X", fit_fct=(lowest_mid_axis, max_along_axis))
        pot_height = ModelParam("pot_height", "Y", fit_fct=(max_along_axis, max_along_axis))
        soil_width = ModelParam("soil_width", "X", fit_fct=(max_normal_axis, lowest_mid_axis))
        soil_height = ModelParam("soil_height", "Y", fit_fct=(lowest_mid_axis, lowest_mid_axis))

        # Set relations between parameters
        pot_height.set_bound(limit_min=soil_width, limit_max=pot_width)
        pot_width.set_bound(limit_max=pot_height)
        soil_width.set_bound(limit_min=pot_height, limit_max=soil_height)
        soil_height.set_bound(limit_max=soil_width)

        # Add all parameters to a list
        self.params = [pot_width, pot_height, soil_width, soil_height]

        # Initialise supporting segments and corresponding point list
        self.segments:list[Segment] = []
        for _ in self.params:
            self.segments.append((Vector((0, 0)), Vector((0, 0))))

        # Initialise model performance
        self.fit_ratio: float = 0

    def __repr__(self) -> str:
        """Output model parameters"""
        model_param = ", ".join(str(param) for param in self.params)
        return f"Model ({model_param})"

    def dist_to_point(self, point:Vector) -> tuple[float, int]:
        """Compute the distance from the point to each segments, return value and segment index"""
        min_dist = float("inf")
        min_index = -1
        for (index, segment) in enumerate(self.segments):
            dist = point_to_segment(point, segment)
            if dist < min_dist:
                min_dist = dist
                min_index = index
        # Return minium distance and matching segment index
        return (min_dist, min_index)

    def init_points_cluster(self) -> list[list[Vector]]:
        """Reset points to segment list"""
        points_cluster:list[list[Vector]] = []
        for _ in self.params:
            points_cluster.append([])
        return points_cluster

    def cluster_points(self, points:list[Vector], dist_thresh:float) -> list[list[Vector]]:
        """Assign points to closest segment is distance is below a threshold"""
        # Init point cluster
        points_cluster = self.init_points_cluster()
        for point in points:
            distance, segment_index = self.dist_to_point(point)
            if distance < dist_thresh:
                # Add point to list corresponding to segment
                points_cluster[segment_index].append(point)
        # Update fitted ratio based on sum of all points assigned to cluster
        self.fit_ratio = sum(len(pts) for pts in points_cluster) / len(points)
        return points_cluster

    def fit_param(self, sampled_points:list[Vector]) -> None:
        """Fit model based on current sampled points"""
        # Fit parameters
        for param in self.params:
            param.set_values(sampled_points)
        # Update supporting segment for all parameters
        for (index, param) in enumerate(self.params):
            self.segments[index] = param.supporting_segment()

    def update_param(self, points_cluster:list[list[Vector]], min_points:int) -> None:
        """Update the parameters based on point cluster for each parameter"""
        for (index, param) in enumerate(self.params):
            # Get relevant coordinate for param and list of points
            axis = param.direction.lower()
            points = points_cluster[index]

            # If enough point, update param value as average of points along axis
            if len(points) >= min_points:
                param.values = param.update_fct(getattr(point, axis) for point in points)

        # Update supporting segment for all parameters
        for (index, param) in enumerate(self.params):
            self.segments[index] = param.supporting_segment()

# RANSAC process to fit model parameters
# TODO: remove deepcopy and dissociate parameters from model
def ransac_fit(model:PotSection,
               points:list[Vector],
               ransac_param:RansacParam) -> PotSection:
    """RANSAC loop, select points from point list, update parameter"""
    # Initialise best model
    best_fit = 0
    best_model = deepcopy(model)

    # Iteratively fit the model to a random sample of points
    for iteration in range(ransac_param.max_iter):
        # Sample points at random
        sampled_points = choices(points, k=ransac_param.nb_sample)
        # Fit new model based on sample points
        sampled_model = deepcopy(model)
        sampled_model.fit_param(sampled_points)
        # Cluster points based on current value of parameters
        points_cluster = sampled_model.cluster_points(points, ransac_param.dist_thresh)

        # Evaluate fitting performance
        # If not points per cluster lower than defined in ransac param, do not update the parameters
        if any(len(cluster) < ransac_param.min_pts_per_line for cluster in points_cluster):
            continue
        # If less outliers than best model, do not update the parameters
        if sampled_model.fit_ratio < best_fit:
            continue

        # Update parameter based on points cluster
#        model.update_param(points_cluster, ransac_param.min_cluster)
        best_fit = sampled_model.fit_ratio
        best_model = sampled_model

        # Evaluate fitting performance
        if best_fit > ransac_param.max_fit:
            print(f"Model fitted in {iteration} iterations")
            return best_model

    # Reach max number of iteration
    print(f"Reaching max {ransac_param.max_iter} iterations")
    return best_model
