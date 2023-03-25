# coding=utf-8
"""Postprocess Toolbox

.. sectionauthor:: Menzenbach, Ruppert
"""

from typing import Tuple, Union, Type
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import interp1d

from pyrit.mesh import Mesh, TriMesh, TetMesh
from pyrit.region import Regions


def smooth(points_on_line: np.ndarray, values: np.ndarray):
    """
    Smoothes unevenly spaced data using spline interpolation and a moving average filter.

    Parameters
    ----------
    points_on_line : np.ndarray
        (N,2) or (N, 3) matrix with coordinates of the evaluated points.
    values : np.ndarray
        (N,) Vector with the value of the component of the field at the evaluated points.

    Returns
    -------
    smoothed_values: np.ndarray
        Vector with the smoothed data.

    Raises
    ------
    ValueError
        ValueError if values does not have dimension (N,) or (N,1)
    ValueError
        ValueError if the numbers of points and values do not match.
    """
    if np.shape(values)[0] != np.size(values):
        raise ValueError("Values must have dimension (N,) or (N,1).")

    if np.shape(values)[0] != np.shape(points_on_line)[0]:
        raise ValueError("The number of values and query points are not the same.")

    x = np.linalg.norm(points_on_line, axis=1)  # Transform points into 1d data

    if np.size(np.unique(np.diff(x))) == 1:  # Points are evenly spaced
        values_smoothed = moving_average_filter(values)
    else:  # Points are not evenly spaced
        # if the number of values is too small for spline interpolation, just return original values
        if np.size(values) < 4:
            return values

        # Create evenly spaced samples using spline interpolation
        spline = interp1d(x, values, kind='cubic')
        x_evenly_sampled = np.linspace(x[0], x[-1], len(x), endpoint=True)
        values_evenly_sampled = spline(x_evenly_sampled)

        # Smooth the evenly sampled data
        y_evenly_sampled = moving_average_filter(values_evenly_sampled)

        # Go back to original query points using spline interpolation
        spline = interp1d(x_evenly_sampled, y_evenly_sampled)
        values_smoothed = spline(x)

    return values_smoothed


def moving_average_filter(values_to_smooth: np.ndarray, window_length: int = 5):
    """
    Smoothes (evenly spaced!) data using a moving average filter with adjusted window size for boundary values \
    (similar to matlab smooth function).

    Parameters
    ----------
    values_to_smooth: np.ndarray
        Evenly spaced data points to be smoothed.
    window_length: int
        window_length of moving average filter.

    Returns
    -------
    smoothed_values: np.ndarray
        Array with the smoothed data.

    Raises
    ------
    ValueError
        ValueError if the window length is an even number.
    ValueError
        ValueError if the data is not one-dimensional, i.e. if the dimension is not (N,) or (N,1)
    """
    # Make sure the window_length is odd.
    if (window_length % 2) == 0:
        raise ValueError("The window_length must be an odd number.")
    if np.shape(values_to_smooth)[0] != np.size(values_to_smooth):
        raise ValueError("Values must have dimension (N,) or (N,1).")

    # General smoothing using a moving average filter with window length=5
    smoothed_values = uniform_filter1d(values_to_smooth, size=window_length)

    # Special treatment for the values at the boundaries
    for idx in range(0, int((window_length - 1) / 2)):
        shortened_window_length = 2 * idx + 1
        smoothed_values[idx] = np.average(values_to_smooth[0:shortened_window_length])
        smoothed_values[-1 - idx] = np.average(values_to_smooth[-shortened_window_length:])

    return smoothed_values


def get_field_on_line(msh: Mesh, start_point: Union[np.ndarray, list, Tuple], end_point: Union[np.ndarray, list, Tuple],
                      num_steps: int, field_values: np.ndarray, smoothing: bool = False,
                      regions: Regions = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns the value of the field in the given direction on the line between the two given points.

    Parameters
    ----------
    msh : Mesh
        The mesh.
    start_point : Union[np.ndarray, list, Tuple]
        Start point of the line which should be evaluated.
    end_point : Union[np.ndarray, list, Tuple]
        End point of the line which should be evaluated.
    num_steps : int
        Number of points on the line which should be evaluated. The points are distributed uniformly on the line.
    field_values : np.ndarray
        Vector with values of one component of the field.
    smoothing : bool
        If set to true the values are smoothend by a Savitzky-Golay filter implemented in scipy:
        :py:func:'scipy.signal.savgol_filter'.
    regions: Regions
        The regions.

    Returns
    -------
    points_on_line : np.ndarray
        (N,2) or (N, 3) matrix with coordinates of the evaluated points.
    values : np.ndarray
        Vector with the value of the component of the field at the evaluated points.

    Raises
    ------
    Exception
        If a given point might not be in the mesh or if another error occurs
        while looking for element/node indices.
    Exception
        If the number of values does not match the number of nodes or elements.
    TypeError
        If the coordinates are 2D but msh is not a 2D TriMesh.
    TypeError
        If the coordinates are 3D but msh is not a 3D TetMesh.
    TypeError
        If smoothing=True but regions is not passed to the function (still None).
    ValueError
        If the coordinates do not have either the dimension 2 or 3.
    """
    start_point_vec = np.asarray(start_point).reshape((1, -1))
    end_point_vec = np.asarray(end_point).reshape((1, -1))
    if len(start_point_vec[0]) == 2 and len(end_point_vec[0]) == 2:
        if not isinstance(msh, TriMesh):  # if not 2D
            raise TypeError("The dimension of the mesh does not match the dimension of the coordinates (2D).")
    elif len(start_point_vec[0]) == 3 and len(end_point_vec[0]) == 3:
        if not isinstance(msh, TetMesh):  # if not 3D
            raise TypeError("The dimension of the mesh does not match the dimension of the coordinates (3D).")
    else:
        raise ValueError("The coordinates of start- and end-point do not have either the dimension 2 or 3.")
    points_on_line = np.concatenate(
        [start_point_vec + i / num_steps * (end_point_vec - start_point_vec) for i in range(num_steps)], axis=0)
    try:
        if len(field_values) == msh.num_node:
            if isinstance(msh, TriMesh):  # if 2D
                # interpolation
                x = msh.node[:, 0]
                y = msh.node[:, 1]
                triangles = msh.elem2node
                triObj = Triangulation(x, y, triangles)
                interpolator = LinearTriInterpolator(triObj, field_values)
                values = interpolator(*points_on_line.T)
            if isinstance(msh, TetMesh):  # if 3D
                raise NotImplementedError("Cannot handle allocation on nodes in 3D yet.")
            elem_idcs = msh.find_elemidx(points_on_line).astype(int)  # for potential smoothing
        elif len(field_values) == msh.num_elem:
            elem_idcs = msh.find_elemidx(points_on_line).astype(int)
            # take only one point if multiple points are located in the same element
            new_indcs = np.argsort(elem_idcs)  # sort to avoid errors when adding unique counts
            sorted_elem_idcs = elem_idcs[new_indcs]
            _, unique_idx, unique_counts = np.unique(sorted_elem_idcs, return_index=True, return_counts=True)
            unique_idx += np.floor(0.5 * unique_counts).astype(int)
            # now the order of the points has to be restored again
            old_idx_order = np.sort(new_indcs[unique_idx])  # [np.where(new_indcs == idx)[0][0] for idx in unique_idx]
            elem_idcs = elem_idcs[old_idx_order]
            points_on_line = points_on_line[old_idx_order, :]
            values = field_values[elem_idcs]
        else:
            raise Exception("The number of values does not match the number of nodes or elements.")
    except Exception as e:
        print(e)
        raise Exception("An Error occured in finding the index of an element/node.\
            Probably one of the points is not in the mesh.") from e
    if smoothing:
        # only smooth over values which belong to the same material to keep jumps in between changes of material
        regions_line = msh.elem2regi[elem_idcs]
        indcs_to_smooth = []
        if regions is None:
            raise TypeError("regions is None. If smoothing is True, regions has to be passed to the function.")
        regions_msh = regions
        last_region = None
        for i, region in enumerate(regions_line):
            if i == 0:
                indcs_to_smooth.append(i)
                last_region = region
            else:
                last_material = regions_msh.get_regi(last_region).mat
                current_material = regions_msh.get_regi(region).mat
                last_region = region
                if last_material != current_material or i == len(regions_line) - 1:
                    values[indcs_to_smooth] = smooth(points_on_line[indcs_to_smooth], values[indcs_to_smooth])
                    values[indcs_to_smooth] = smooth(points_on_line[indcs_to_smooth], values[indcs_to_smooth])

                    indcs_to_smooth = []  # reset indc list
                    indcs_to_smooth.append(i)  # add current index (current index belongs to new material)
                else:
                    indcs_to_smooth.append(i)
    return points_on_line, values


def plot_field_on_line(msh: Mesh, start_point: Union[np.ndarray, list, Tuple],
                       end_point: Union[np.ndarray, list, Tuple], field_values: np.ndarray,
                       lower_bound: float = 0, upper_bound: float = 1, x_label: str = 'x axis',
                       y_label: str = 'y axis', title: str = '', smoothing: bool = False,
                       regions=None, fig: Figure = None, ax: Type[Axes] = None, **kwargs) -> Tuple[Figure, Type[Axes]]:
    """
    Plots the value of the field in the given direction on the line between the two given points.

    The function :py:func:`get_field_on_line` gets called with the following defaults:
    num_steps = 200
    window_length = 3
    polyorder = 1

    Parameters
    ----------
    msh : Mesh
        The mesh.
    start_point : Union[np.ndarray, list, Tuple]
        Start point of the line which should be evaluated.
    end_point : Union[np.ndarray, list, Tuple]
        End point of the line which should be evaluated.
    field_values : np.ndarray
        Vector with values of one component of the field.
    lower_bound : float
        Lower bound for first axis. Deafault is 0.
    upper_bound : float
        Upper bound for first axis. Deafault is 1.
    x_label : str, optional
        Label of the x axis. Default is 'x axis'
    y_label : str, optional
        Label of the y axis. Default is 'y axis'
    title : str, optional
        Title of the plot. Default is ''
    smoothing : bool, optional
        If set to true the values are smoothend by a savgol-filter. Default is None
    regions : Regions, optional
        A regions object. Default is None
    fig : Figure, optional
        A figure object. Default is None
    ax : Axes, optional
        An axes object. Default is None


    Returns
    -------
    fig: Figure
        An object of matplotlib.figure.Figure.
    ax: Axes
        An object of matplotlib.axes.Axes

    Raises
    ------
    Exception
        If a given point might not be in the mesh or if another error occurs
        while looking for element/node indices.
    Exception
        If the number of values does not match the number of nodes or elements.
    TypeError
        If the coordinates are 2D but msh is not a 2D TriMesh.
    TypeError
        If the coordinates are 3D but msh is not a 3D TetMesh.
    ValueError
        If the coordinates do not have either the dimension 2 or 3.
    ValueError
        If lower_bound is >= upper_bound.
    """
    start_point = np.asarray(start_point).reshape((1, -1))
    end_point = np.asarray(end_point).reshape((1, -1))

    points_on_line, values = get_field_on_line(msh, start_point, end_point, 200, field_values, smoothing, regions)
    # determine values of the first axis
    total_length = np.linalg.norm(start_point - end_point)
    rel_distances_points = np.linalg.norm(points_on_line - start_point, axis=1) / total_length
    if lower_bound >= upper_bound:
        raise ValueError("Lower_bound is >= upper_bound but should be lower.")
    first_axis = lower_bound + (upper_bound - lower_bound) * rel_distances_points

    if fig is None or ax is None:  # either both or none should be None
        fig = plt.figure()
        ax = fig.add_subplot(111)

    ax.plot(first_axis, values, **kwargs)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    plt.title(title)
    return fig, ax
