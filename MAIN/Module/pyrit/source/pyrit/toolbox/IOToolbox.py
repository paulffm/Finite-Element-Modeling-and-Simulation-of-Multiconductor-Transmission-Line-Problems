# coding=utf-8
"""Input-Output-Toolbox

.. sectionauthor:: Bundschuh
"""

from typing import Union, List, Any, Iterable, Dict, Tuple, TYPE_CHECKING
from pathlib import Path
from contextlib import contextmanager
import pickle
import datetime as dt
import io
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

from pyrit import get_logger

logger = get_logger(__name__)

if TYPE_CHECKING:
    from pyrit.mesh import TriMesh


def save(obj: Any, path: Union[str, 'Path'], ignore_attributes: List[str] = None) -> 'Path':
    """Saves the given object.

    Parameters
    ----------
    obj : Any
        Any object that should be saved in a file.
    path : Union[str, Path]
        The path. If a string is given, it is converted to a Path object. The file ending can but must not be given.
        In any case it will be a *.pkl* file. See `pickle doc <https://docs.python.org/3/library/pickle.html>`_ for
        more information.
        If the path contains a non-existing folder, this folder will be created.
    ignore_attributes: List[str], optional
        A list of attributes that are not included in the saved file. By default, all attributes are included.

    Returns
    -------
    path : Path
        The path that was used for opening the file.
    """
    path = Path(path).with_suffix('.pkl')

    if not path.parent.exists():
        logger.info("Creating path: %s.", str(path.parent))
        path.parent.mkdir()

    tmp_dict = {}
    if ignore_attributes:
        for attribute in ignore_attributes:
            try:
                tmp_dict[attribute] = obj.__getattribute__(attribute)
                delattr(obj, attribute)
            except AttributeError:
                logger.warning("The object does not contain the attribute '%s'. It is ignored.", attribute)

    logger.info("Saving to %s.", str(path))
    with open(path, 'wb') as file:
        pickle.dump(obj, file, pickle.HIGHEST_PROTOCOL)
    logger.info("Done saving.")

    if ignore_attributes:
        for attribute in ignore_attributes:
            try:
                obj.__setattr__(attribute, tmp_dict[attribute])
            except KeyError:
                logger.info("The key '%s' does not exist and the object attribute cannot be set.", attribute)

    return path


def load(path: Union[str, Path]) -> Any:
    """Load an object from a file.

    Parameters
    ----------
    path : Union[str, Path]
        The path. If a string is given, it is converted to a Path object. The file ending can but must not be given.
        In any case it will be a *.pkl* file. See `pickle doc <https://docs.python.org/3/library/pickle.html>`_ for
        more information.

    Returns
    -------
    obj : Any
        The loaded object.
    """
    path = Path(path).with_suffix('.pkl')

    logger.info("Start reading problem from %s.", str(path))
    with open(path, 'rb') as file:
        obj = pickle.load(file)
    logger.info("Done reading problem.")
    return obj


class ToLaTeX:
    """Export data for plotting in LaTeX using TikZ and pgfplots."""

    def __init__(self, **kwargs):
        """Constructor.

        Parameters
        ----------
        kwargs :
            Default keyword arguments for subsequent methods. If a keyword argument is not given in a method call, but
            was given in the constructor, it is used in the method. So you can set a default for multiple methods.
        """
        self.kwargs = kwargs

    def _parse_kwargs(self, **kwargs):
        """Parse the kwargs.

        If a kwarg is not given in a function, it is taken from the object (if given).
        """
        for key, val in self.kwargs.items():
            if key not in kwargs:
                kwargs[key] = val
        return kwargs

    def comment(self, msg: str, **kwargs) -> str:
        """Convert a string to a comment.

        Splits the `msh` at line breaks and starts every line with the `comment_string` of the `kwargs`.

        Parameters
        ----------
        msg : str
            The message to process.
        kwargs :
            From the keyword arguments, the key 'comment_string' is used.

        Returns
        -------
        out : str
            The formatted message.
        """
        kwargs = self._parse_kwargs(**kwargs)

        comment_string = kwargs.setdefault('comment_string', '# ')

        msg = msg.strip()
        out = ''
        for s in msg.split('\n'):
            out = out + comment_string + s + '\n'
        return out

    @staticmethod
    def check_parent(path: Path):
        """Creates the parent directory of the path, if it does not exist.

        Parameters
        ----------
        path : Path
        """
        if not path.parent.exists():
            logger.info("Creating path: %s.", str(path.parent))
            path.parent.mkdir()

    @contextmanager
    def open(self, path: Path, message: Union[str, Iterable[str]] = None, **kwargs):
        """Context manager for opening a file and adding basic information.

        Parameters
        ----------
        path : Path
            The path.
        message : Union[str, Iterable[str]], optional
            The message to add at the top of the file.
        kwargs :
            Additional keyword arguments. The key 'file_header' is used.
        """
        kwargs = self._parse_kwargs(**kwargs)

        file_header = kwargs.setdefault('file_header', True)
        logger.info("Saving to %s.", str(path))
        with open(path, 'wt') as file:  # pylint: disable=unspecified-encoding
            if isinstance(file_header, bool):
                if file_header:
                    file.write(self.comment('Created: ' + dt.datetime.now().strftime("%B %d, %Y; %H:%M:%S"), **kwargs))
            else:
                file.write(file_header)
            if message is not None:
                if isinstance(message, str):
                    message = [message, ]
                for m in message:
                    file.write(self.comment(m, **kwargs))
            yield file
        logger.info("Done saving.")

    def save_dat(self, data: Union[pd.DataFrame, Dict, Iterable[np.ndarray], np.ndarray], path: Union[str, Path],
                 message: Union[str, Iterable[str]] = None, **kwargs) -> Path:
        r"""Save data to a .dat file.

        Save data to .dat file. This file format can be used to plot the data with pgfplots in LaTeX.
        Besides the data, you can give a number of messages that are put at the top of the file as a comment.

        Parameters
        ----------
        data : Union[pd.DataFrame, Dict, Iterable[np.ndarray], np.ndarray]
            The data to save. If ``data`` is a DataFrame, nothing is changed. If not, it is converted to a DataFrame.
        path : Union[str, Path]
            The path where to save the file. If a string is given, it is converted to a Path object. The file ending can
            but must not be given.
        message : Union[str, Iterable[str]]
            A message put as a comment at the top of the file.
        kwargs : Any
            Additional arguments. With 'comment_string' you can determine the string that is used to initialize a
            commented line (used for messages).
            The remaining arguments are passed to :py:meth:`pandas.DataFrame.to_csv`

        Returns
        -------
        path : Path
            The path where the file was saved.

        Examples
        --------
        Here a small example how to use the generated .dat-file to plot in LaTeX.

        Suppose you saved following data:

        >>> time = np.arange(10)
        >>> data = {'time': time, 'func': time**2}
        >>> save_dat(data, 'data')

        Then, the most simple way to plot this data in LaTeX is by:

        .. code-block:: latex

            \begin{tikzpicture}
                \begin{axis}[xlabel={time $t$},ylabel={$f(t)$}]
                    \addplot+[] table[x=time, y=func] {data.dat};
                \end{axis}
            \end{tikzpicture}

        Note that the packages `tikzpicture` and `pgfplots` have to be loaded. This plot can then be extended as
        desired. For more information about this, see https://www.ctan.org/pkg/pgf
        """
        kwargs = self._parse_kwargs(**kwargs)

        path = Path(path).with_suffix('.dat')

        if isinstance(data, pd.DataFrame):
            df = data
        else:
            df = pd.DataFrame(data)

        sep = kwargs.pop('sep', ' ')
        line_terminator = kwargs.pop('line_terminator', '\n')
        index = kwargs.pop('index', False)

        with self.open(path, message, **kwargs) as file:
            kwargs.pop('comment_string', None)
            kwargs.pop('file_header', None)
            file.write(df.to_csv(sep=sep, line_terminator=line_terminator, index=index, **kwargs))

        return path

    def load_dat(self, path: Union[str, Path], **kwargs) -> pd.DataFrame:
        """Load from .dat-file.

        Parameters
        ----------
        path : Union[str, Path]
            The path to load from.
        kwargs :
            Additional arguments. With `comment_string` you can determine which lines should be ignored (when they begin
            with that string). The remaining arguments are passed to :py:meth:`pandas.DataFrame.load_csv`.

        Returns
        -------
        data : pd.DataFrame
            The loaded data
        """
        kwargs = self._parse_kwargs(**kwargs)

        path = Path(path).with_suffix('.dat')

        sep = kwargs.pop('sep', ' ')
        comment_string = kwargs.pop('comment_string', '#')
        logger.info("Start reading .dat-file from %s.", str(path))
        with open(path, 'rt') as file:  # pylint: disable=unspecified-encoding
            lines = file.readlines()
            for k, line in enumerate(lines):
                if not line.startswith(comment_string):
                    break
            data_lines = lines[k:]  # pylint: disable=undefined-loop-variable
            data = pd.read_csv(io.StringIO(''.join(data_lines)), sep=sep, **kwargs)
        logger.info("Done reading.")

        return data

    def export_trimesh(self, mesh: 'TriMesh', path: Union[str, Path], message: Union[str, Iterable[str]] = None,
                       one_file: bool = False, **kwargs) -> Union[Path, Tuple[Path]]:
        r"""Save a trimesh to a dat file.

        Parameters
        ----------
        mesh : TriMesh
            The mesh.
        path : Path
            The path.
        message : Union[str, Iterable[str]], optional
            The message to add at the top of the file. If an Iterable is given. Each element is a separate line.
            Default is None
        one_file : bool, optional
            If True, one file is generated. If False, two files are generated. Two files use less memory.
        kwargs :
            Additional keyword arguments. 'sep' is the separator between two entries, 'line_terminator' is the string
            that ends a line. If `one_file` is False, 'nodes_file_appendix' is a string that is added to the file name
            of the file with node information and 'elements_file_appendix' to the file name of the file with element
            information. The rest is passed to subsequent method calls (:py:meth:`comment` and :py:meth:`open`).

        Returns
        -------
        path: Union[Path, Tuple[Path]]
            If `one_file` is True, the path of this file. If `one_file` is False, the path of the elements and the path
            of the nodes.

        Notes
        -----
        The pgfplotslibrary *patchplots* can be loaded for more options, but is not necessary.
        The LaTeX code for one file is:

        .. code-block:: latex

            \begin{tikzpicture}
                \begin{axis}[axis equal image, colormap={whitecm}{color=(white) color=(white)}]
                    \addplot[line width=0.1pt, faceted color=black] [patch] table {mesh.dat};
                \end{axis}
            \end{tikzpicture}

        and for two files:

        .. code-block:: latex

            \begin{tikzpicture}%two files
                \begin{axis}[axis equal image, colormap={whitecm}{color=(white) color=(white)}]
                    \addplot[line width=0.1pt,faceted color=black, patch table={mesh_elements.dat}] [patch] table
                    {mesh_nodes.dat};
                \end{axis}
            \end{tikzpicture}
        """
        kwargs = self._parse_kwargs(**kwargs)

        path = Path(path).with_suffix('.dat')

        sep = kwargs.pop('sep', ' ')
        line_terminator = kwargs.pop('line_terminator', '\n')

        if one_file:
            with self.open(path, message, **kwargs) as file:
                for nodes in mesh.elem2node:
                    for node in nodes:
                        n = mesh.node[node]
                        file.write(f"{n[0]}" + sep + f"{n[1]}" + line_terminator)
            return path
        # else

        nodes_file_appendix = kwargs.pop('nodes_file_appendix', '_nodes')
        elements_file_appendix = kwargs.pop('elements_file_appendix', '_elements')

        path_node = path.with_name(path.stem + nodes_file_appendix).with_suffix(path.suffix)
        path_element = path.with_name(path.stem + elements_file_appendix).with_suffix(path.suffix)

        with self.open(path_node, message, **kwargs) as file:
            file.write(self.comment(f'This is the node file. '
                                    f'The corresponding element file is \'{path_element.name}\'.', **kwargs))
            for node in mesh.node:
                file.write(f"{node[0]}" + sep + f"{node[1]}" + line_terminator)

        with self.open(path_element, message, **kwargs) as file:
            file.write(self.comment(f'This is the element file. '
                                    f'The corresponding node file is \'{path_node.name}\'.', **kwargs))
            for element in mesh.elem2node:
                file.write(f"{element[0]}" + sep + f"{element[1]}" + sep + f"{element[2]}" + line_terminator)

        return path_node, path_element

    def export_field(self, mesh: 'TriMesh', field: np.ndarray, path: Union[str, Path],
                     message: Union[str, Iterable[str]] = None, one_file: bool = False, **kwargs) \
            -> Union[Path, Tuple[Path]]:
        r"""Export a field to a .dat file.

        A field that is allocated on the elements or the nodes of the mesh is exported.

        Parameters
        ----------
        mesh : TriMesh
            The mesh
        field : np.ndarray
            The field. Must be a one dimensional array. Its size must be the number of elements or the number of nodes
            of the mesh.
        path : Union[str, Path]
            A Path or a string representing the path. The suffix will be ignored.
        message : Union[str, Iterable[str]], optional
            A message or a number of messages placed at the top of the written files. Default is None
        one_file : bool
            If True, only one file is generated. If False, two files are generated. This uses less memory.
        kwargs :
            Additional keyword arguments. 'sep' is the separator between two entries, 'line_terminator' is the string
            that ends a line. If `one_file` is False, 'nodes_file_appendix' is a string that is added to the file name
            of the file with node information and 'elements_file_appendix' to the file name of the file with element
            information. The rest is passed to subsequent method calls (:py:meth:`comment` and :py:meth:`open`).

        Returns
        -------
        path: Union[Path, Tuple[Path]]
            If `one_file` is True, the path of this file. If `one_file` is False, the path of the elements and the path
            of the nodes.

        Notes
        -----
        In the following code examples, only the field is plotted. If you also want to plot the mesh, add the options
        `line width=0.1pt,faceted color=black` in the `\addplot`-command. The pgfplotslibrary *patchplots* can be loaded
        for more options, but is not necessary.

        The LaTeX code for one file (nodal data and element wise data) is:

        .. code-block:: latex

            \begin{tikzpicture}
                \begin{axis}[axis equal image,colormap name=viridis, colorbar]
                    \addplot [patch,shader=interp]
                    table [point meta=\thisrow{v}] {file.dat};
                \end{axis}
            \end{tikzpicture}

        For nodal data and two files:

        .. code-block:: latex

            \begin{tikzpicture}
                \begin{axis}[axis equal image, colormap name=viridis, colorbar]
                    \addplot [line width=0pt,faceted color=none,patch,patch table={file_elements.dat}]
                    table [point meta=\thisrow{v}] {file_nodes.dat};
                \end{axis}
            \end{tikzpicture}

        For element wise data and two files:

        .. code-block:: latex

            \begin{tikzpicture}
                \begin{axis}[axis equal image, colormap name=viridis, colorbar]
                    \addplot [line width=0pt,faceted color=none, patch,patch table with point meta={file_elements.dat}]
                    table {file_nodes.dat};
                \end{axis}
            \end{tikzpicture}
        """
        kwargs = self._parse_kwargs(**kwargs)

        path = Path(path).with_suffix('.dat')

        flag_nodal: bool
        if field.shape[0] == mesh.num_node:
            flag_nodal = True
        elif field.shape[0] == mesh.num_elem:
            flag_nodal = False
        else:
            raise ValueError(f"field has length {field.shape[0]} but should be {mesh.num_node}, the number of nodes.")

        sep = kwargs.pop('sep', ' ')
        line_terminator = kwargs.pop('line_terminator', '\n')

        if one_file and flag_nodal:
            with self.open(path, message, **kwargs) as file:
                file.write("r" + sep + "z" + sep + "v" + line_terminator)
                for nodes in mesh.elem2node:
                    for node in nodes:
                        n = mesh.node[node]
                        file.write(f"{n[0]}" + sep + f"{n[1]}" + sep + f"{field[node]}" + line_terminator)

            return path

        if one_file and not flag_nodal:
            with self.open(path, message, **kwargs) as file:
                file.write("r" + sep + "z" + sep + "v" + line_terminator)
                for element, nodes in enumerate(mesh.elem2node):
                    for node in nodes:
                        n = mesh.node[node]
                        file.write(f"{n[0]}" + sep + f"{n[1]}" + sep + f"{field[element]}" + line_terminator)

            return path

        # one_file is False

        nodes_file_appendix = kwargs.pop('nodes_file_appendix', '_nodes')
        elements_file_appendix = kwargs.pop('elements_file_appendix', '_elements')

        path_node = path.with_name(path.stem + nodes_file_appendix).with_suffix(path.suffix)
        path_element = path.with_name(path.stem + elements_file_appendix).with_suffix(path.suffix)

        if flag_nodal:
            # file with elements to nodes
            with self.open(path_element, message, **kwargs) as file:
                file.write(self.comment(f'This is the element file. '
                                        f'The corresponding node file is \'{path_node.name}\'.', **kwargs))
                for element in mesh.elem2node:
                    file.write(f"{element[0]}" + sep + f"{element[1]}" + sep + f"{element[2]}" + line_terminator)

            # file with nodes and values
            with self.open(path_node, message, **kwargs) as file:
                file.write(
                    self.comment(f'This is the node file. '
                                 f'The corresponding element file is \'{path_element.name}\'.', **kwargs))
                file.write("r" + sep + "z" + sep + "v" + line_terminator)
                for k, node in enumerate(mesh.node):
                    file.write(f"{node[0]}" + sep + f"{node[1]}" + sep + f"{field[k]}" + line_terminator)
        else:
            # file with elements and values
            with self.open(path_element, message, **kwargs) as file:
                file.write(self.comment(f'This is the element file. '
                                        f'The corresponding node file is \'{path_node.name}\'.', **kwargs))
                for k, element in enumerate(mesh.elem2node):
                    file.write(f"{element[0]}" + sep + f"{element[1]}" + sep + f"{element[2]}" + sep +  # noqa:W504
                               f"{field[k]}" + sep + line_terminator)

            # file with nodes
            with self.open(path_node, message, **kwargs) as file:
                file.write(
                    self.comment(f'This is the node file. '
                                 f'The corresponding element file is \'{path_element.name}\'.', **kwargs))
                for node in mesh.node:
                    file.write(f"{node[0]}" + sep + f"{node[1]}" + line_terminator)

        return path_element, path_node

    def export_nodes(self, mesh: 'TriMesh', node_indices: Iterable[int], path: Union[str, Path],
                     message: Union[str, Iterable[str]] = None, **kwargs) -> Path:
        """Export nodes of the mesh to LaTeX.

        Parameters
        ----------
        mesh : TriMesh
            The mesh.
        node_indices : Iterable[int]
            Iterable ove node indices.
        path : Union[str, Path]
            A Path or a string representing the path. The suffix will be ignored.
        message : Union[str, Iterable[str]]
            A message or a number of messages placed at the top of the written files. Default is None
        kwargs :
            Additional keyword arguments. 'tikz_cs' is the coordinate system, the points are declared in, 'options'
            is a string of options that is added to the draw command, 'command' is the tikz-command for a line and
            'radius' is the radius of a node. The rest is passed to :py:meth:`open`.

        Returns
        -------
        path : Path
            The path of the saved file.
        """
        kwargs = self._parse_kwargs(**kwargs)

        path = Path(path).with_suffix('.tex')

        if 'comment_string' not in kwargs:
            kwargs['comment_string'] = '% '

        tikz_cs = kwargs.pop('tikz_cs', 'axis cs:')
        command = kwargs.pop('command', r"\fill")
        options = kwargs.pop('options', "black")
        radius = kwargs.pop('radius', '1pt')

        with self.open(path, message, **kwargs) as file:
            for index in node_indices:
                file.write(f"{command}[{options}] ({tikz_cs}{mesh.node[index][0]},{mesh.node[index][1]}) "
                           f"circle ({radius});\n")

        return path

    def export_edges(self, mesh: 'TriMesh', edge_indices: Iterable[int], path: Union[str, Path],
                     message: Union[str, Iterable[str]] = None, **kwargs) -> Path:
        """Export edges of the mesh to LaTeX

        Parameters
        ----------
        mesh : TriMesh
            The mesh.
        edge_indices : Iterable[int]
            Iterable over edge indices.
        path : Union[str, Path]
            A Path or a string representing the path. The suffix will be ignored.
        message : Union[str, Iterable[str]], optional
            A message or a number of messages placed at the top of the written files. Default is None
        kwargs :
            Additional keyword arguments. 'tikz_cs' is the coordinate system, the points are declared in, and 'options'
            is a string of options that is added to the draw command. The rest is passed to :py:meth:`open`.

        Returns
        -------
        path : Path
            The path of the saved file.
        """
        kwargs = self._parse_kwargs(**kwargs)

        path = Path(path).with_suffix('.tex')
        tikz_cs = kwargs.pop('tikz_cs', 'axis cs:')
        options = kwargs.pop('options', None)

        if 'comment_string' not in kwargs:
            kwargs['comment_string'] = '% '

        if options is not None:
            option_str = f"[{options}]"
        else:
            option_str = ''

        with self.open(path, message, **kwargs) as file:
            for index in edge_indices:
                node1, node2 = mesh.edge2node[index]
                file.write(r"\draw" + f"{option_str} ({tikz_cs}{mesh.node[node1][0]},{mesh.node[node1][1]}) -- "
                                      f"({tikz_cs}{mesh.node[node2][0]},{mesh.node[node2][1]});\n")

        return path

    def export_elements(self, mesh: 'TriMesh', element_indices: Iterable[int], path: Union[str, Path],
                        message: Union[str, Iterable[str]] = None, **kwargs) -> Path:
        """Export elements of the mesh to LaTeX

        Parameters
        ----------
        mesh : TriMesh
            The mesh.
        element_indices : Iterable[int]
            Iterable over the element indices.
        path : Path
            A Path or a string representing the path. The suffix will be ignored.
        message : Union[str, Iterable[str]], optional
            A message or a number of messages placed at the top of the written files. Default is None
        kwargs :
            Additional keyword arguments. 'tikz_cs' is the coordinate system, the points are declared in, and 'options'
            is a string of options that is added to the draw command. The rest is passed to :py:meth:`open`.

        Returns
        -------
        path : Path
            The path of the saved file.
        """
        kwargs = self._parse_kwargs(**kwargs)

        path = Path(path).with_suffix('.tex')

        if 'comment_string' not in kwargs:
            kwargs['comment_string'] = '% '

        tikz_cs = kwargs.pop('tikz_cs', 'axis cs:')
        options = kwargs.pop('options', None)

        if options is not None:
            option_str = f"[{options}]"
        else:
            option_str = ''

        with self.open(path, message, **kwargs) as file:
            for index in element_indices:
                node1, node2, node3 = mesh.elem2node[index]
                file.write(r"\draw " + f"{option_str}" +  # noqa:W504
                           f"({tikz_cs}{mesh.node[node1][0]},{mesh.node[node1][1]}) -- "
                           f"({tikz_cs}{mesh.node[node2][0]},{mesh.node[node2][1]}) -- "
                           f"({tikz_cs}{mesh.node[node3][0]},{mesh.node[node3][1]}) -- cycle;\n")

        return path

    def export_region_edges(self, mesh: 'TriMesh', path: Union[str, Path], message: Union[str, Iterable[str]] = None,
                            **kwargs) -> Path:
        """Export all edges that divide two regions or are at the outer boundary.

        Use :py:meth:`export_edges` to export the edges.

        Parameters
        ----------
        mesh : TriMesh
            The mesh.
        path : Union[str, Path]
            A Path or a string representing the path. The suffix will be ignored.
        message : Union[str, Iterable[str]], optional
            A message or a number of messages placed at the top of the written files. Default is None
        kwargs :
            Additional keyword arguments. All passed to :py:meth:`export_edges`.

        Returns
        -------
        path : Path
            The path of the saved file.
        """
        kwargs = self._parse_kwargs(**kwargs)

        if 'comment_string' not in kwargs:
            kwargs['comment_string'] = '% '

        edges = []
        for edge in range(mesh.num_edge):
            elems = mesh.edge2elem[edge]
            if mesh.elem2regi[elems[0]] != mesh.elem2regi[elems[1]] or elems[1] == -1:
                edges.append(edge)

        return self.export_edges(mesh, edges, path, message, **kwargs)

    def export_contour(self, mesh: 'TriMesh', field: np.ndarray, path: Union[str, Path],
                       message: Union[str, Iterable[str]] = None, **kwargs):
        r"""Export the contour lines of a field.

        Parameters
        ----------
        mesh : TriMesh
            A mesh object
        field : np.ndarray
            A field
        path : Union[str, Path]
            The path to save to the file to
        message : Union[str, Iterable[str]]
            A message that it put at the top of the file. If an iterable, each item is put on a separate line. Each line
            is started with the comment string.
        kwargs :
            Other options.

        Notes
        -----
        Instruction to display the result in LaTeX:

        .. code-block:: latex

            \begin{tikzpicture}
                \begin{axis}[axis equal image, unbounded coords=jump]
                    \addplot table {data.dat};
                \end{axis}
            \end{tikzpicture}

        Returns
        -------
        path : Path
            The path of the saved file.
        """
        kwargs = self._parse_kwargs(**kwargs)
        num_lines = kwargs.pop('num_Lines', 40)
        path = Path(path).with_suffix('.dat')

        fig, ax, _ = mesh.plot_equilines(field, num_lines)
        plt.close(fig)
        lines = []
        for col in ax.collections:
            try:
                segments = col.get_paths()
                for segment in segments:
                    lines.append(segment.vertices)
            except AttributeError:
                print("Attribute Error")

        # lines = [col.get_segments()[0] for col in ax.collections]

        sep = kwargs.pop('sep', ' ')
        line_terminator = kwargs.pop('line_terminator', '\n')
        with self.open(path, message, **kwargs) as file:
            # file.write("r" + sep + "z" + sep + "v" + line_terminator)
            for line in lines:
                x, y = line[:, 0], line[:, 1]
                num_points = len(x)
                for k in range(num_points):
                    file.write(str(x[k]) + sep + str(y[k]) + line_terminator)
                file.write(str(x[0]) + sep + str(y[0]) + line_terminator)
                file.write('0' + sep + 'nan' + line_terminator)

        return path
