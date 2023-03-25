# coding=utf-8
"""Field Plotter

.. sectionauthor:: Bundschuh
"""

from typing import TYPE_CHECKING, Tuple, Type, Iterable, List
from abc import ABC
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.artist import Artist
from matplotlib import animation

from pyrit import get_logger

if TYPE_CHECKING:
    from pyrit.mesh import TriMesh
    from matplotlib.animation import ArtistAnimation

logger = get_logger(__name__)

__all__ = ['FieldPlotter', 'StaticFieldPlotter', 'HarmonicFieldPlotter', 'TransientFieldPlotter']


class FieldPlotter(ABC):
    """A class to generalize the plotting of fields on triangular meshes"""

    def __init__(self, tri_mesh: 'TriMesh'):
        """The Constructor.

        Parameters
        ----------
        tri_mesh : TriMesh
            A mesh object that will be used for the plotting.
        """
        self.mesh = tri_mesh

    def plot_scalar_field(self, field: np.ndarray, **kwargs) -> Tuple[Figure, Axes, Artist]:
        """Plot a scalar field.

        Parameters
        ----------
        field : np.ndarray
            A field to plot.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`.

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        artist : Artist
            The artist
        """
        return self.mesh.plot_scalar_field(field, **kwargs)

    def plot_vector_field(self, field: np.ndarray, plot_type: str, **kwargs) -> Tuple[Figure, Axes, Artist]:
        """Plot a vector field.

        Parameters
        ----------
        field : np.ndarray
            The field to plot.
        plot_type : str
            The plot_type of the plot. One of 'arrows', 'abs', 'stream'.
        kwargs :
            Passed to the plot method.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        artist : Artist
            The artist
        """
        if plot_type == 'arrows':
            return self.mesh.plot_vector_field(field, **kwargs)
        if plot_type == 'abs':
            return self.mesh.plot_scalar_field(np.linalg.norm(field, axis=1), **kwargs)
        if plot_type == 'stream':
            return self.mesh.plot_streamlines(field, **kwargs)
        logger.warning("The given plot_type '%s' is not known. The plot_type 'arrows' is used instead.", plot_type)
        return self.mesh.plot_vector_field(field, **kwargs)

    def plot_equilines(self, field: np.ndarray, **kwargs) -> Tuple[Figure, Axes, Artist]:
        """Plot equilines of the field.

        Parameters
        ----------
        field : np.ndarray
            A field to plot.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_equilines`.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        artist : Artist
            The artist
        """
        return self.mesh.plot_equilines(field, **kwargs)


class StaticFieldPlotter(FieldPlotter):
    """A field plotter for static fields."""


class HarmonicFieldPlotter(FieldPlotter):
    """A field plotter for harmonic fields.

    The fields given fields are considered as phasors, i.e. should be complex valued.
    """

    def __init__(self, mesh: 'TriMesh', frequency: float):
        """The constructor.

        Parameters
        ----------
        mesh : TriMesh
            A mesh object that will be used for the plotting.
        frequency : float
            The frequency of the harmonic fields.
        """
        super().__init__(mesh)
        self.frequency = frequency

    @staticmethod
    def get_field_at_phase(frequency: float, field: np.ndarray, phase: float = 0, time: float = None):
        r"""Evaluate the complex valued field at a certain phase.

        Each entry in field is a complex number, representing a phasor. Return the field at a certain phase, or for a
        certain time point. Let :math:`\varphi` be the phase and :math:`\underline{g}` be a complex value. Then the
        formula is

        .. math::

            g = \Re\left\{\underline{g} e^{\jmath\varphi}\right\}\,.

        If the time :math:`t` is given, the phase is determined by :math:`\varphi = 2\pi f t`, where :math:`f` is the
        frequency.

        Parameters
        ----------
        frequency : float
            The frequency of the field.
        field : np.ndarray
            The field.
        phase : float, optional
            The phase. Default is 0
        time : float, optional
            The time. Default is None

        Returns
        -------
        field_at_phase : np.ndarray
            Field at the specified phase or time.
        """
        if time is not None:
            logger.info("Using the variable 'time' instead of 'phase'.")
            phase = np.mod(2 * np.pi * frequency * time, 2 * np.pi)
        field_at_phase = np.real(field * np.exp(1j * phase))
        return field_at_phase, phase

    @staticmethod
    def update_title(phase: float, **kwargs):
        """Updates the title. Adds the phase information.

        Searches the keyword 'title'. If it exists, append the phase information. If not, create it with phase
        information.

        Parameters
        ----------
        phase : float
            The phase
        kwargs :
            Keyword arguments.

        Returns
        -------
        kwargs :
            Keyword arguments with ensured 'title' that contains the phase.
        """
        old_title = kwargs.get('title', '')
        if old_title == '':
            title = f'Phase {phase}'
        else:
            title = old_title + f' (at phase {phase})'
        kwargs['title'] = title

        return kwargs

    def process_inputs(self, field: np.ndarray, **kwargs):
        """Process inputs for harmonic filed plotting.

        Determines the field at phase and updates the title.

        Parameters
        ----------
        field : np.ndarray
            The field
        kwargs :
            Keyword arguments. Keys 'phase' and 'time' are processed here and removed. The rest is passed to the
            plotting methods.

        Returns
        -------
        field_at_phase : np.ndarray
            The field at the given phase
        kwargs :
            The processed kwargs.
        """
        phase = kwargs.pop('phase', 0)
        time = kwargs.pop('time', None)

        # Determine the field at the right phase
        field_at_phase, phase = self.get_field_at_phase(self.frequency, field, phase, time)

        # Update the title
        kwargs = self.update_title(phase, **kwargs)

        return field_at_phase, kwargs

    def plot_scalar_field(self, field, **kwargs):
        field_at_phase, kwargs = self.process_inputs(field, **kwargs)

        return super().plot_scalar_field(field_at_phase, **kwargs)

    def plot_vector_field(self, field: np.ndarray, plot_type: str, **kwargs) -> Tuple[Figure, Type[Axes]]:
        field_at_phase, kwargs = self.process_inputs(field, **kwargs)

        return super().plot_vector_field(field_at_phase, plot_type, **kwargs)

    def plot_equilines(self, field, **kwargs):
        field_at_phase, kwargs = self.process_inputs(field, **kwargs)

        return super().plot_equilines(field_at_phase, **kwargs)


class TransientFieldPlotter(FieldPlotter):
    """A field plotter for transient solutions."""

    @staticmethod
    def update_title(time: float, **kwargs):
        """Update the title in the kwargs

        Parameters
        ----------
        time :
        kwargs :

        Returns
        -------
        kwargs :
            The updated kwargs.
        """
        if time is None:
            return kwargs
        old_title = kwargs.get('title', '')
        if old_title == '':
            title = f'Time {time:.4g}'
        else:
            title = old_title + f' (at time {time:.4g})'
        kwargs['title'] = title

        return kwargs

    def plot_scalar_field(self, field, time=None, **kwargs):
        kwargs = self.update_title(time, **kwargs)
        return super().plot_scalar_field(field, **kwargs)

    def plot_vector_field(self, field: np.ndarray, plot_type: str, time=None, **kwargs):
        kwargs = self.update_title(time, **kwargs)
        return super().plot_vector_field(field, plot_type, **kwargs)

    def plot_equilines(self, field, time=None, **kwargs):
        kwargs = self.update_title(time, **kwargs)
        return super().plot_equilines(field, **kwargs)

    @staticmethod
    def plot_over_time_scalar(time_steps: np.ndarray, values: np.ndarray, **kwargs) -> Tuple[Figure, Axes, Artist]:
        """Plot a number of scalar value over time.

        Parameters
        ----------
        time_steps : np.ndarray
            An array of the time steps.
        values : np.ndarray
            An array with the values to plot over time.
        kwargs :
            Arguments passed to the plot function.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        line : Artist
            The plotted line.
        """
        fig, ax = plt.subplots()
        ax: Axes
        ax.set_title(kwargs.pop('title', ''))
        ax.set_xlabel(kwargs.pop('xlabel', 'time'))
        ax.set_ylabel(kwargs.pop('ylabel', ''))
        line = ax.plot(time_steps, values, **kwargs)
        return fig, ax, line

    @staticmethod
    def plot_over_time_vector(time_steps: np.ndarray, values: np.ndarray, **kwargs) -> \
            Tuple[Figure, Axes, List[Artist]]:
        """Plot a number of vectors over time.

        Parameters
        ----------
        time_steps : np.ndarray
            An array of the time steps.
        values : np.ndarray
            An array with the values to plot over time. The first dimension of the array must be equal to the length of
            the time steps.
        kwargs :
            Arguments passed to the plot function.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        lines : List[Artist]
            A list of the plotted lines.
        """
        fig = plt.figure()
        ax: Axes = fig.add_subplot(111, projection='3d')
        len_data = values.shape[1]
        x = kwargs.pop('x', np.arange(len_data))
        ax.set_title(kwargs.pop('title', ''))
        ax.set_xlabel(kwargs.pop('xlabel', ''))
        ax.set_ylabel(kwargs.pop('ylabel', 'time'))
        lines = []
        for k, ts in enumerate(time_steps):
            line = ax.plot(x, ts * np.ones(len_data), values[k], **kwargs)
            lines.append(line)
        return fig, ax, lines

    def plot_over_time(self, time_steps: np.ndarray, values: np.ndarray, **kwargs):
        """Plot the specified values over time.

        If value is a 2D array, the method :py:meth:`plot_over_time_vector` is used. If it is a 1D array the method
        :py:meth:`plot_over_time_scalar` is ised.

        Parameters
        ----------
        time_steps : np.ndarray
            Array with time steps
        values : np.ndarray
            Array with values to plot
        kwargs :
            Arguments passed to the called method

        Returns
        -------
        Return of the called method
        """
        if len(values.shape) > 1:
            return self.plot_over_time_vector(time_steps, values, **kwargs)
        return self.plot_over_time_scalar(time_steps, values, **kwargs)

    def animate_scalar_field(self, fields: Iterable[np.ndarray], x_label: str = 'x axis', y_label: str = 'y axis',
                             title: str = '', **kwargs) -> Tuple[Figure, Axes, 'ArtistAnimation']:
        """Animate a scalar field.

        Parameters
        ----------
        fields : Iterable[np.ndarray]
            An iterable over the fields at the time instances
        x_label : str, optional
            The x label. Default is 'x axis'.
        y_label : str, optional
            The y label. Default is 'y axis'.
        title : str, optional
            The title. Default is ''
        kwargs :
            Additional settings

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        animation : ArtistAnimation
            The animation object
        """
        animation_kwargs = kwargs.pop('animation_kwargs', {})
        animation_kwargs.setdefault('interval', 50)
        animation_kwargs.setdefault('blit', True)
        animation_kwargs.setdefault('repeat_delay', 1000)

        fig = kwargs.pop('fig', None)
        ax = kwargs.pop('ax', None)
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
        ax: Axes
        artists = []
        for k, field in enumerate(fields):
            if k == 0:
                self.plot_scalar_field(field, fig=fig, ax=ax, just_plot=True)
            _, _, artist = self.plot_scalar_field(field, fig=fig, ax=ax, just_plot=True)
            artists.append([artist])

        ani = animation.ArtistAnimation(fig, artists, **animation_kwargs)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        plt.colorbar(artists[-1][0], shrink=0.5, aspect=5)
        plt.show()

        return fig, ax, ani

    def animate_vector_field(self, fields: Iterable[np.ndarray], plot_type: str, x_label: str = 'x axis',
                             y_label: str = 'y axis', title: str = '', **kwargs) -> \
            Tuple[Figure, Axes, 'ArtistAnimation']:
        """Animate a vector field.

        Parameters
        ----------
        fields : Iterable[np.ndarray]
            An iterable over the fields at the time instances
        plot_type : str
            The plot_type of the plot. See :py:meth:`plot_vector_field`
        x_label : str, optional
            The x label. Default is 'x axis'.
        y_label : str, optional
            The y label. Default is 'y axis'.
        title : str, optional
            The title. Default is ''
        kwargs :
            Additional settings

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        animation : ArtistAnimation
            The animation object
        """
        animation_kwargs = kwargs.pop('animation_kwargs', {})
        animation_kwargs.setdefault('interval', 50)
        animation_kwargs.setdefault('blit', True)
        animation_kwargs.setdefault('repeat_delay', 1000)

        fig = kwargs.pop('fig', None)
        ax = kwargs.pop('ax', None)
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)

        artists = []
        for k, field in enumerate(fields):
            if k == 0:
                self.plot_vector_field(field, plot_type, fig=fig, ax=ax, just_plot=True)
            _, _, artist = self.plot_vector_field(field, plot_type, fig=fig, ax=ax, just_plot=True)
            artists.append([artist])

        ani = animation.ArtistAnimation(fig, artists, **animation_kwargs)

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        plt.colorbar(artists[-1][0], shrink=0.5, aspect=5)

        plt.show()

        return fig, ax, ani

    def animate_equilines(self, fields, num_lines: int = 20, x_label: str = 'x axis', y_label: str = 'y axis',
                          title: str = '', **kwargs):
        """Animate equilines.

        Parameters
        ----------
        fields : Iterable[np.ndarray]
            An iterable over the fields at the time instances
        num_lines : int, optional
            The number of equilines. Default is 20.
        x_label : str, optional
            The x label. Default is 'x axis'.
        y_label : str, optional
            The y label. Default is 'y axis'.
        title : str, optional
            The title. Default is ''
        kwargs :
            Additional settings

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        animation : ArtistAnimation
            The animation object
        """
        animation_kwargs = kwargs.pop('animation_kwargs', {})
        animation_kwargs.setdefault('interval', 50)
        animation_kwargs.setdefault('blit', True)
        animation_kwargs.setdefault('repeat_delay', 1000)

        fig = kwargs.pop('fig', None)
        ax = kwargs.pop('ax', None)
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)

        # fields = [field for field in fields]
        artists = []
        for k, field in enumerate(fields):
            if k == 0:
                self.plot_equilines(field, num_lines, fig=fig, ax=ax, just_plot=True)
            _, _, artist = self.plot_equilines(field, num_lines, fig=fig, ax=ax, just_plot=True)
            artists.append([artist])

        ani = animation.ArtistAnimation(fig, artists, **animation_kwargs)

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        plt.colorbar(artists[-1][0], shrink=0.5, aspect=5)

        plt.show()

        return fig, ax, ani
