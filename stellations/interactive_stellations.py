""" Interactive tool to view Hamiltonian Edge Colourings for planar graphs.
    This tool barely works. It crudely uses interactive matplotlib to enable the
    real-time interaction of graphs, and thus this program has several issues.
    Namely, if Qt5agg does not work on your machine, you may need to tweak which
    renderer matplotlib uses, and that might require fiddling with other parts of
    the program.
    Author: Max Croucher
    Date: 13/10/24
"""


from sys import argv
import subprocess
from copy import deepcopy
import os
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Button
import networkx as nx
from graphs import Graph
from read_pl import read_graphs, PLWriter, write_graph

matplotlib.use("Qt5agg")

FRAMERATE = 20

# Default Filenames for the required binaries
if os.name == 'posix':
    CYCLE_GENERATOR = Path(__file__).absolute().parent.joinpath(Path("bin", "all_ham_cycles"))
    PATH_GENERATOR = Path(__file__).absolute().parent.joinpath(Path("bin", "all_longest_paths"))
elif os.name == 'nx':
    CYCLE_GENERATOR = Path(__file__).absolute().parent.joinpath(Path("bin", "all_ham_cycles.exe"))
    PATH_GENERATOR = Path(__file__).absolute().parent.joinpath(Path("bin", "all_longest_paths.exe"))
else:
    raise NotImplementedError(f"Operating System {os.name} is not supported")
TEMP = Path(os.getcwd()).joinpath(Path("temp"))
TEMP.mkdir(parents=True, exist_ok=True)

FRAME_COORDS = [0.1, 0.15, 0.8, 0.8]

BUTTON_WIDTH = 0.19
BUTTON_HEIGHT = 0.03
BUTTON_XPOS = [0.1, 0.3033333, 0.5066666, 0.71]
BUTTON_YPOS = [0.1, 0.05]
RECTANGLE_PADDING = 0.0075

BUTTON_COLOUR = {
    "default": (0.85,0.85,0.85),
    "valid": (0.8,1,0.8),
    "invalid": (1,0.7,0.7),
    "active": (0.65,1,0.65),
    "inactive": (1,0.65,0.65)
}


class ButtonWrapper():
    """Wrapper for a matplotlib Button widget, to better handle position, states and events"""
    def __init__(self, pos, text, states=None, active=None, is_valid=True, col=None, hov_col=None):
        self.pltfig = None
        self.axis = None
        self.button = None
        self.pos = pos
        self.text = text
        if states is None:
            self.states = (False, True)
        else:
            self.states = states
        self._stateref = 0

        if col is None:
            self.colours = ('active', 'inactive')
        else:
            self.colours = col

        if hov_col is None:
            self.hover_colours = ('valid', 'invalid')
        else:
            self.hover_colours = hov_col

        if active is None:
            self.active_states = [0]
        else:
            self.active_states = active
        self.is_valid = is_valid


    def update_button(self):
        """Update display of button object"""
        self.button.color = BUTTON_COLOUR[
            self.colours[0] if self._stateref in self.active_states else self.colours[-1]
        ]
        self.button.hovercolor = BUTTON_COLOUR[
            self.hover_colours[0] if self.is_valid else self.hover_colours[-1]
        ]
        self.button.label.set_text(self.text)


    def event(self, _):
        """event handler to be invoked when clicking the button"""
        self._stateref = (self._stateref + 1) % len(self.states)


    def create_button(self, pltfig):
        """Initialise a button object in the correct figure"""
        self.pltfig = pltfig
        self.axis = self.pltfig.add_axes([
            BUTTON_XPOS[self.pos[0]],
            BUTTON_YPOS[self.pos[1]],
            BUTTON_WIDTH, BUTTON_HEIGHT
        ])
        self.button = Button(self.axis, self.text)
        self.update_button()
        self.button.on_clicked(self.event)


    @property
    def state(self):
        """Get the current button state"""
        return self.states[self._stateref]


    def set_state(self, new_state):
        """Update the current button state"""
        if new_state not in self.states:
            raise ValueError("Invalid State")
        self._stateref = self.states.index(new_state)


    def set_text(self, text):
        """Set the button's display text"""
        self.text = text


def get_longest_path_length():
    """Find the length of the longest path in the current graph"""
    with open(TEMP.joinpath("temp.paths"), encoding="utf-8") as path_obj:
        while '}\n' != path_obj.readline():
            pass
        return len(path_obj.readline().split(' - ')) - 1


def is_ccw(vec_a, vec_b, vec_c):
    """Return true if a triangle abc is counter-clockwise."""
    vec_p = vec_b - vec_a
    vec_q = vec_c - vec_a
    area = vec_p[0] * vec_q[1] - vec_q[0] * vec_p[1]
    return area > 0


def mean_center(polygon):
    """Find the mean-center of a triangle"""
    return sum(v for v in polygon) / len(polygon)


def perimeter(coords, face):
    """Return the perimeter of a face"""
    return sum(
        np.linalg.norm(coords[face[i]]-coords[face[(i+1)%len(face)]])
     for i in range(len(face)))


EVENT_STATES = {
    'CLICK_EVENT': False,
    'CLICK_STATE': False,
    'VALID_MOUSE': False,
    'HOME_BUTTON': False,
    'IN_FRAME': False,
    'CLICK_POS': np.array([0, 0]),
    'MOUSE_POS': np.array([0, 0])
}

BUTTONS = {
    "colouring": ButtonWrapper((0, 0),'Toggle Gradient Colouring'),
    "back": ButtonWrapper((2, 0),'Back (0)',col=['default']),
    "forward": ButtonWrapper((3, 0),'Forward (0)',col=['default']),
    "redraw": ButtonWrapper((0, 1),'Redraw Graph',col=['default']),
    "freqency": ButtonWrapper((1, 0),'Toggle Frequencies',states=['Off','Count','Freq'],active=[1,2]),
    "stellate": ButtonWrapper((1, 1),'Stellate',states=['Off','Clicked','On'],active=[1,2]),
    "delete": ButtonWrapper((2, 1),'Delete',states=['Off','Clicked','On'],active=[1,2]),
    "permute": ButtonWrapper((3, 1),'Change Outer Face',states=['Off','Clicked','On'],active=[1,2])
}
linked_BUTTONS = [['stellate', 'delete', 'permute']]
BUTTONS['stellate'].set_state('On')


def get_faces(adj):
    """Return a set of the faces of a graph"""
    faces = set()
    for v in adj:
        for u in adj[v]:
            w = adj[u][(adj[u].index(v) + 1) % len (adj[u])]
            current_face = [v, u]
            while w != v:
                current_face.append(w)
                next_vertex = adj[w][(adj[w].index(u) + 1) % len (adj[w])]
                u = w
                w = next_vertex
            min_index = current_face.index(min(current_face))
            faces.add(tuple(current_face[min_index:] + current_face[:min_index]))
    return faces


def get_coloured_edges(path_file_name, count_cycles=False, get_frequencies=False, path_mode="H"):
    """Use a pre-generated path file to find the set of edges that lie on some Hamiltonian Cycle"""
    if path_mode == "N":
        output = [{}]
        if count_cycles:
            output.append(0)
        if get_frequencies:
            output.append({})
        yield output
    else:
        with open(path_file_name, encoding='utf-8') as path_file:
            num_cycles = 0
            freqs = dict()
            coloured_edges = dict()
            def next_line():
                _curr = path_file.readline().strip()
                if _curr == '':
                    raise StopIteration
                return _curr
            try:
                current_line = next_line()
                while True:
                    found_cycles = set()
                    num_cycles = 0
                    if current_line != "graph = {":
                        raise ValueError(f"Error. Expected 'graph = {{'. Got {current_line}")
                    while current_line != "}":
                        current_line = next_line()
                    current_line = next_line()
                    coloured_edges = dict()
                    while current_line != "graph = {":
                        vertices = current_line.split(' - ')
                        if get_frequencies:
                            min_v_index = vertices.index(min(vertices))
                            if path_mode == "L":
                                signature = tuple(vertices)
                            else:
                                signature = tuple(vertices[min_v_index:] + vertices[:min_v_index])
                            if signature in found_cycles:
                                current_line = next_line()
                                continue
                            rev_vert = list(reversed(vertices))
                            min_v_index = rev_vert.index(min(rev_vert))
                            if path_mode == "L":
                                signature = tuple(rev_vert)
                            else:
                                signature = tuple(rev_vert[min_v_index:] + rev_vert[:min_v_index])
                            if signature in found_cycles:
                                current_line = next_line()
                                continue
                            found_cycles.add(signature)
                        num_cycles += 1
                        for i in range(path_mode=="L", len(vertices)):
                            e = tuple(sorted([int(vertices[i-1])+1, int(vertices[i])+1]))
                            freqs[e] = freqs.get(e, 0) + 1
                            coloured_edges[e] = 'r'
                        current_line = next_line()

                outdata = [coloured_edges]
                if count_cycles:
                    outdata.append(num_cycles)
                if get_frequencies:
                    outdata.append(freqs)
                yield tuple(outdata)
            except StopIteration:
                outdata = [coloured_edges]
                if count_cycles:
                    outdata.append(num_cycles)
                if get_frequencies:
                    outdata.append(freqs)
                yield tuple(outdata)
            except ValueError:
                return


def get_common_edges(path_file_name, path_mode="L"):
    """Use a pre-generated path file to find the set of edges that lie on every Hamiltonian Cycle"""
    if path_mode == "N":
        yield {}
    else:
        with open(path_file_name, encoding='utf-8') as path_file:
            common_edges = None
            def next_line():
                _curr = path_file.readline().strip()
                if _curr == '':
                    raise StopIteration
                return _curr
            try:
                current_line = next_line()
                while True:
                    if current_line != "graph = {":
                        raise ValueError(f"Error. Expected 'graph = {{'. Got {current_line}")
                    while current_line != "}":
                        current_line = next_line()
                    current_line = next_line()
                    common_edges = None
                    while current_line != "graph = {":
                        vertices = current_line.split(' - ')
                        edges = {tuple(sorted(
                            [int(vertices[i-1])+1, int(vertices[i])+1]
                        ))for i in range(path_mode=="L", len(vertices))}
                        if common_edges is None:
                            common_edges = edges
                        else:
                            common_edges &= edges
                        current_line = next_line()
                    yield set() if common_edges is None else common_edges
            except StopIteration:
                yield set() if common_edges is None else common_edges
            except ValueError:
                return


def get_endpoints(path_file_name):
    """Use a pre-generated path file to find the set of vertices that are endpoints of some path"""
    with open(path_file_name, encoding='utf-8') as path_file:
        all_endpoints = None
        def next_line():
            _curr = path_file.readline().strip()
            if _curr == '':
                raise StopIteration
            return _curr
        try:
            current_line = next_line()
            while True:
                if current_line != "graph = {":
                    raise ValueError(f"Error. Expected 'graph = {{'. Got {current_line}")
                while current_line != "}":
                    current_line = next_line()
                current_line = next_line()
                all_endpoints = None
                while current_line != "graph = {":
                    vertices = current_line.split(' - ')
                    endpoints = set([int(vertices[0])+1, int(vertices[-1])+1])
                    if all_endpoints is None:
                        all_endpoints = endpoints
                    else:
                        all_endpoints |= endpoints
                    current_line = next_line()
                yield set() if all_endpoints is None else all_endpoints
        except StopIteration:
            yield set() if all_endpoints is None else all_endpoints
        except ValueError:
            return


def closest_point(coords, point):
    """Return the closest point from a set of points"""
    return min((np.linalg.norm(v-point), k) for k, v in coords.items())


def distance_to_edge(edge, point):
    """returns the distance between a line segment and a point"""
    vec_u, vec_v = edge
    vec_a = point - vec_u
    vec_b = vec_v - vec_u
    factor = (np.dot(vec_a,vec_b) / np.dot(vec_b,vec_b))
    if factor <= 0:
        return np.linalg.norm(vec_u-point)
    if factor >= 1:
        return np.linalg.norm(vec_v-point)
    proj = vec_u + vec_b * factor
    return np.linalg.norm(proj-point)


def closest_edge(coords, adj, point):
    """Returns the closest edge from a set of points"""
    return min(
        min(
            (distance_to_edge((coords[u],coords[v]), point), tuple(sorted((u, v)))) for u in adj[v])
        for v in adj
    )


def is_intersecting(edge_1, edge_2):
    """Determines whether two line segments intersect"""
    u1, v1 = edge_1
    u2, v2 = edge_2
    return is_ccw(u1, v1, u2) != is_ccw(u1, v1, v2) and is_ccw(u2, v2, u1) != is_ccw(u2, v2, v1)


def vec_angle(vector):
    """Return the angle between vector and [0,1]"""
    sign = (-1 if vector[1]<0 else 1)
    dotprod = np.dot(vector,np.array([1, 0]) * (1 / np.linalg.norm(vector)))
    return np.rad2deg(np.arccos(dotprod))*sign


def onclick(event):
    """Handle the left mouse button being clicked"""
    EVENT_STATES['CLICK_POS'] = (
        event.xdata if event.xdata is not None else EVENT_STATES['CLICK_POS'][0],
        event.ydata if event.ydata is not None else EVENT_STATES['CLICK_POS'][1]
    )
    EVENT_STATES['CLICK_EVENT'] = True
    EVENT_STATES['CLICK_STATE'] = True


def offclick(event):
    """Handle the left mouse button being unclicked"""
    if event.xdata is not None and event.ydata is not None:
        EVENT_STATES['CLICK_POS'] = np.array([event.xdata, event.ydata])
    EVENT_STATES['CLICK_STATE'] = False


def mouse_move(event):
    """Retrieve current mouse position"""
    EVENT_STATES['VALID_MOUSE'] = all([
        event.xdata is not None,
        event.ydata is not None,
        EVENT_STATES['IN_FRAME']
    ])
    if EVENT_STATES['VALID_MOUSE']:
        EVENT_STATES['MOUSE_POS'] = np.array([event.xdata, event.ydata])


def home_event(_):
    """Handle the home button being clicked"""
    EVENT_STATES['HOME_BUTTON'] = True


def axis_enter(event):
    """Handle the mouse entering a frame"""
    axis_position = event.inaxes.get_position()
    EVENT_STATES['IN_FRAME'] = all([
        axis_position.x0 == FRAME_COORDS[0],
        axis_position.y0 == FRAME_COORDS[1]
    ])


def axis_exit(_):
    """Handle the mouse leaving a frame"""
    EVENT_STATES['IN_FRAME'] = False


def plural(integer):
    """return an 's' if integer is not 1"""
    return '' if integer == 1 else 's'


def scale_unit(coords):
    """Scale a set of points on a 2d plane to lie between (-1, -1) and (1, 1)"""
    min_pos = np.min(list(coords.values()), axis=0)
    max_pos = np.max(list(coords.values()), axis=0)
    midpoint = np.mean([min_pos,max_pos], axis=0)
    spread = np.abs(midpoint-min_pos)
    return {v: (pos-midpoint)/spread for v, pos in coords.items()}


def get_contained_face(coords, faces, point):
    """Find the (possibly non-convex) face that surrounds a point"""
    max_pos = np.max(list(coords.values()), axis=0)
    trace_line = (point, max_pos)
    for face in faces:
        num_intersections = sum(
            is_intersecting(trace_line, (coords[face[i]], coords[face[(i+1)%len(face)]]))
        for i in range(len(face)))
        if num_intersections % 2:
            return face
    return faces[-1]


def generate_edge_data(adj, path_mode="L"):
    """Get edge colouring data from an adjacency list.
       Involves storing the graph to a temporary file"""
    if path_mode != "N":
        with PLWriter(TEMP.joinpath("temp.pl")) as pl_obj:
            write_graph(pl_obj, adj)
        with open(TEMP.joinpath("temp.paths"), 'w', encoding="utf-8") as path_obj:
            process = subprocess.Popen(
                [(PATH_GENERATOR if path_mode=="L" else CYCLE_GENERATOR),str(TEMP.joinpath("temp.pl"))],
                stdout=path_obj
            )
            process.wait()
    coloured_edges, num_cycles, e_freqs = next(get_coloured_edges(
        TEMP.joinpath("temp.paths"),
        count_cycles=True,
        get_frequencies=True,
        path_mode=path_mode
    ))
    common_edges = next(
        get_common_edges(TEMP.joinpath("temp.paths"), path_mode=path_mode)
    )
    print(f"{len(coloured_edges)} red")
    print(f"{len(common_edges)} blue")
    for edge in common_edges:
        coloured_edges[edge] = 'cyan'
    if path_mode != "L":
        endpoints=set()
    else:
        endpoints = next(get_endpoints(TEMP.joinpath("temp.paths")))
    return coloured_edges, num_cycles, e_freqs, endpoints


def draw_from_coords(adj, v_pos, colours=None, e_freqs=None, num_cycles=None, do_endpoints=None):
    """Generate Matplotlib objects for edges and vertices of a graph"""
    if colours is None:
        colours = dict()
    if do_endpoints is None:
        do_endpoints = set()
    v_obj = plt.scatter(
        *zip(*list(v_pos.values())),
        150,
        color=['g' if v in do_endpoints else 'k' for v in v_pos],
        marker='o',
        zorder=10
    )
    t_objs = dict()
    for v, (x, y) in v_pos.items():
        t_objs[v] = plt.text(
            x,y,v,
            horizontalalignment='center',
            verticalalignment='center',
            color='w',
            fontsize='medium',
            fontweight='bold',
            zorder=20
        )
    e_objs = dict()
    e_text_objs = dict()
    for v, neighbours in adj.items():
        for u in neighbours:
            edge_colour = colours.get((u, v), colours.get((v, u), 'k'))
            if isinstance(edge_colour, tuple):
                colour_tup = edge_colour
            else:
                colour_tup = matplotlib.colors.to_rgb(edge_colour)
            inverse_colour = tuple(1-c for c in colour_tup)
            if tuple(sorted([v, u])) not in e_objs:
                if BUTTONS['colouring'].state:
                    if e_freqs.get(tuple(sorted([v, u])), 0) != num_cycles:
                        edge_colour = (e_freqs.get(tuple(sorted([v, u])), 0)/num_cycles, 0, 0)
                e_objs[tuple(sorted([v, u]))] = plt.plot(
                    *zip(v_pos[v], v_pos[u]),
                    color=edge_colour, linewidth=3, zorder=0
                )[0]
            if BUTTONS['freqency'].state != 'Off':
                count = e_freqs.get(tuple(sorted([v, u])), 0)
                if BUTTONS['freqency'].state == 'Freq':
                    display_val = f"{count/num_cycles:.3f}"
                else:
                    display_val = count
                if tuple(sorted([v, u])) not in e_text_objs and 0 < count < num_cycles:
                    mean_x = (v_pos[u][0] + v_pos[v][0])/2
                    mean_y = (v_pos[u][1] + v_pos[v][1])/2
                    e_text_objs[tuple(sorted([v, u]))] = plt.text(
                        mean_x,mean_y,display_val,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color=inverse_colour,
                        backgroundcolor=colour_tup,
                        fontsize='medium',
                        zorder=20
                    )
    return v_obj, t_objs, e_objs, e_text_objs


def draw_interactive(n, output=None, path_mode='H'):
    """Draw an interactive matplotlib window"""
    fig = plt.figure(figsize=(10, 10))

    # widget button handlers
    for button in BUTTONS.values():
        button.create_button(fig)

    #set up axes
    ax = fig.add_axes(FRAME_COORDS)
    ax.set_xticks([])
    ax.set_yticks([])
    g = nx.PlanarEmbedding()
    g.set_data(n)
    plt.xlim((-1.2, 1.2))
    plt.ylim((-1.2, 1.2))

    #UX rectangle
    ux_rect = matplotlib.patches.Rectangle(
        (BUTTON_XPOS[1]-RECTANGLE_PADDING,
        BUTTON_YPOS[1]-RECTANGLE_PADDING),
        #(0.1, 0.15),
        BUTTON_XPOS[3]+BUTTON_WIDTH-BUTTON_XPOS[1]+2*RECTANGLE_PADDING,
        BUTTON_YPOS[1]+BUTTON_HEIGHT-BUTTON_YPOS[1]+2*RECTANGLE_PADDING,
        #0.8,0.8,
        clip_on=False, transform=fig.transFigure, fill=None, ec='k')
    ax.add_patch(ux_rect)

    # Get data about current graph
    faces = list(get_faces(n))
    e_cols, num_cycles, e_freqs, endpoints = generate_edge_data(n, path_mode=path_mode)
    initial_pos = {k: np.array(v) for k, v in nx.planar_layout(g).items()}
    pos = scale_unit(initial_pos)

    # Initialise the output IO if needed
    save_message = ""
    output_index = 0
    if output is not None:
        pl_worker = PLWriter(output)
        pl_object = pl_worker.__enter__()
        write_graph(pl_object, n)
        save_message = f"Most recent graph was saved to {output.name} at index {output_index}"
    if path_mode == 'L':
        path_name = f"Longest Path{plural(num_cycles)} with length {get_longest_path_length()}"
    elif path_mode == "H":
        path_name = f"Hamiltonian Cycle{plural(num_cycles)}"
    else:
        path_name = ''
    save_text = plt.text(
        -1.2,-1.6,
        '' if path_mode == 'N' else f"Current graph has {num_cycles} {path_name}. {save_message}",
        fontsize=12
    )

    #Reorder faces so the outer face is the last one
    outer_face = max(faces, key=lambda x: perimeter(pos, x))
    faces.remove(outer_face)
    faces.append(outer_face)

    #Initialise graph history
    history = [(
        deepcopy(n),
        deepcopy(pos),
        deepcopy(faces),
        deepcopy(e_cols),
        num_cycles,
        deepcopy(endpoints),
        deepcopy(e_freqs)
    )]
    history_index = 0
    v_obj, t_objs, e_objs, e_text_objs = draw_from_coords(
        n,
        pos,
        colours=e_cols,
        e_freqs=e_freqs,
        num_cycles=num_cycles,
        do_endpoints=endpoints
    )

    # Activate event management
    matplotlib_events = {
        'button_manager': plt.get_current_fig_manager(),
        'mouse': fig.canvas.mpl_connect('motion_notify_event', mouse_move),
        'click': fig.canvas.mpl_connect('button_press_event', onclick),
        'unclick': fig.canvas.mpl_connect('button_release_event', offclick),
        'enter_axis': fig.canvas.mpl_connect('axes_enter_event', axis_enter),
        'exit_axis': fig.canvas.mpl_connect('axes_leave_event', axis_exit)
    }
    matplotlib_events['button_manager'].toolbar.actions()[0].triggered.connect(home_event)
    mouse_state = 'F'

    # Initialise interactive highlighting
    plt.title(f"Face: {faces[-1]}")
    face_points = [tuple(pos[faces[-1][i%len(faces[-1])]]) for i in range(len(faces[-1])+1)]
    face_poly = mpatches.Polygon(face_points, fill=True, color='g', alpha=0.3)
    vertex_circle = mpatches.Circle((0, 0), radius=0, color='g', alpha=0.3)
    edge_ellipse = mpatches.Ellipse((0, 0), 0, 0.05, color='g', alpha=0.3)
    ax.add_patch(face_poly)
    ax.add_patch(vertex_circle)
    ax.add_patch(edge_ellipse)

    plt.ion() # Begin interactive loop
    last_active_state = 'stellate'
    while plt.fignum_exists(fig.number):
        for link in linked_BUTTONS:
            for b in link:
                if BUTTONS[b].state == 'Clicked':
                    for x in link:
                        BUTTONS[x].set_state('Off')
                    BUTTONS[b].set_state('On')
                    last_active_state = b
            if all(BUTTONS[x].state == 'Off' for x in link):
                BUTTONS[last_active_state].set_state('On')
        face_poly.set_color('b' if BUTTONS['permute'].state == "On" else 'g')

        # Get mouse info
        face = get_contained_face(pos, faces, EVENT_STATES['MOUSE_POS'])
        face_points = [tuple(pos[face[i%len(face)]]) for i in range(len(face)+1)]
        face_poly.xy = face_points
        vertex_distance, closest_p = closest_point(pos, EVENT_STATES['MOUSE_POS'])
        edge_distance, closest_e = closest_edge(pos, n, EVENT_STATES['MOUSE_POS'])
        vertex_circle.center = pos[closest_p]

        #Collect mouse event and determine which relevant graph feature is closest
        if not EVENT_STATES['CLICK_STATE']:
            if vertex_distance < 0.03 and BUTTONS['permute'].state != "On":
                mouse_state = 'V'
            elif edge_distance < 0.03 and BUTTONS['delete'].state == "On":
                mouse_state = 'E'
            elif BUTTONS['delete'].state == "Off":
                mouse_state = 'F'
            else:
                mouse_state = 'E'

        # Mouse is hovering over a face
        if mouse_state == 'F' and EVENT_STATES['VALID_MOUSE']:
            face_poly.set_alpha(0.3)
            plt.title(f"Selected Face: {face}")
        else:
            face_poly.set_alpha(0)
        if mouse_state == 'E':
            u, v = sorted([pos[closest_e[0]], pos[closest_e[1]]], key=tuple)
            d = u - v
            edge_ellipse.set_center(tuple((u + v) * 0.5))
            edge_ellipse.set_width((np.linalg.norm(u-v)) + 0.05)
            edge_ellipse.set_angle(vec_angle(d))
            edge_ellipse.set_alpha(0.3)
        else:
            edge_ellipse.set_alpha(0)

        # Mouse is hovering over a vertex
        if mouse_state == 'V' and not EVENT_STATES['CLICK_STATE']:
            vertex_circle.radius = 0.05
            plt.title(f"Selected Vertex: {closest_p}")
        else:
            vertex_circle.radius = 0

        # Mouse is clicking and dragging a vertex
        if all([
            mouse_state == 'V',
            EVENT_STATES['CLICK_STATE'],
            BUTTONS['delete'].state == "Off"
        ]):
            plt.title(f"Moving Vertex {closest_p}")
            pos[closest_p] = EVENT_STATES['MOUSE_POS']

        if mouse_state == "E":
            plt.title(f"Selected Edge: {closest_e}")

        # Mouse is not in frame
        if not EVENT_STATES['VALID_MOUSE']:
            plt.title("Mouse out of bounds")

        # Selected object cannot be deleted
        if BUTTONS['delete'].state == 'On':
            if any([
                mouse_state == "E" and any([len(n[closest_e[0]]) <= 2, len(n[closest_e[1]]) <= 2]),
                mouse_state == "V" and any([len(n[u]) <= 2 for u in n[closest_p]])
            ]):
                vertex_circle.set_color('r')
                edge_ellipse.set_color('r')
            else:
                vertex_circle.set_color('g')
                edge_ellipse.set_color('g')

        # Mouse has just been clicked
        if EVENT_STATES['CLICK_EVENT'] and EVENT_STATES['VALID_MOUSE']:
            if BUTTONS['permute'].state == "On" and mouse_state == "F": # Change outer Face
                new_n = dict()
                # Reorder adjacency data:
                # (1) make first vertex in adjacency a vertex of the face
                new_n[face[0]] = deepcopy(n[face[0]])
                # (2) Fill out rest of adjacency
                for v in n:
                    new_n[v] = deepcopy(n[v])
                # (3) Ensure first edge is a boundary edge of the face
                target_index = new_n[face[0]].index(face[1])
                if new_n[face[0]][(target_index-1) % len(new_n[face[0]])] in face:
                    target_index = (target_index-1) % len(new_n[face[0]])
                new_n[face[0]] = new_n[face[0]][target_index:] + new_n[face[0]][:target_index]

                # Get new position data
                g = nx.PlanarEmbedding()
                g.set_data(new_n)
                initial_new_pos = {k: np.array(v) for k, v in nx.planar_layout(g).items()}
                new_pos = scale_unit(initial_new_pos)

                # Recompute faces
                new_faces = list(get_faces(new_n))
                outer_face = max(new_faces, key=lambda x: perimeter(new_pos, x))
                new_faces.remove(outer_face)
                new_faces.append(outer_face)

                print(f"Changed outer face to {face}")
                print(new_n)

                # Update history
                history = history[:history_index+1]
                history.append((
                    deepcopy(new_n),
                    deepcopy(new_pos),
                    deepcopy(new_faces),
                    deepcopy(e_cols),
                    num_cycles,
                    deepcopy(endpoints),
                    deepcopy(e_freqs)
                ))
                history_index += 1
                n,pos,faces,e_cols,num_cycles,endpoints,e_freqs = history[history_index]

            elif BUTTONS['stellate'].state == "On" and mouse_state == "F": # Stellate and add face
                g = Graph()
                g.import_dict(n)
                new_vertex = g.stellate(face)
                new_n = g.get_adjacency()
                c = mean_center([pos[v] for v in face])
                new_pos = deepcopy(pos)
                new_pos[new_vertex] = c

                # Update title
                plt.title("Processing ...")
                plt.draw()
                plt.pause(0.05)

                # Slow edge colouring function
                new_e_cols, new_num_cycles, new_e_freqs, endpoints = generate_edge_data(new_n, path_mode=path_mode)

                # If the stellated face was the outer face
                if face == outer_face:
                    g = nx.PlanarEmbedding()
                    g.set_data(new_n)
                    initial_new_pos = {k: np.array(v) for k, v in nx.planar_layout(g).items()}
                    new_pos = scale_unit(initial_new_pos)

                # Get new graph data
                new_faces = list(get_faces(new_n))
                outer_face = max(new_faces, key=lambda x: perimeter(new_pos, x))
                new_faces.remove(outer_face)
                new_faces.append(outer_face)

                print(f"Stellated {face}")
                print(new_n)

                # Update graph history
                history = history[:history_index+1]
                history.append((
                    deepcopy(new_n),
                    deepcopy(new_pos),
                    deepcopy(new_faces),
                    deepcopy(new_e_cols),
                    new_num_cycles,
                    deepcopy(endpoints),
                    deepcopy(new_e_freqs)
                ))
                history_index+=1
                n, pos, faces, e_cols, num_cycles, endpoints, e_freqs = history[history_index]

                # Export graph
                if output is not None:
                    write_graph(pl_object, n)
                    output_index += 1

            elif BUTTONS['delete'].state == "On":
                new_n = deepcopy(n)
                new_pos = deepcopy(pos)
                valid = True
                if mouse_state == "E":
                    if len(new_n[closest_e[0]]) <= 2 or len(new_n[closest_e[1]]) <= 2:
                        print("Cannot isolate a vertex or introduce bridges!")
                        valid = False
                    else:
                        new_n[closest_e[0]].remove(closest_e[1])
                        new_n[closest_e[1]].remove(closest_e[0])
                        print(f"Deleted Edge {closest_e}")
                if mouse_state == "V":
                    if any([len(new_n[u]) <= 2 for u in new_n[closest_p]]):
                        print("Cannot isolate a vertex or introduce bridges!")
                        valid = False
                    else:
                        for v in new_n[closest_p]:
                            new_n[v].remove(closest_p)
                        del new_n[closest_p]
                        del new_pos[closest_p]
                        max_v = max(new_n)
                        if closest_p < max_v: # Rename vertices to preserve order
                            new_n[closest_p] = new_n[max_v]
                            del new_n[max_v]
                            for u in new_n:
                                if max_v in new_n[u]:
                                    new_n[u][new_n[u].index(max_v)] = closest_p
                            new_pos[closest_p] = new_pos[max_v]
                            del new_pos[max_v]
                        print(f"Deleted Vertex {closest_p}. Renamed vertex {max_v} to {closest_p}")
                if valid:
                    # Update title
                    plt.title("Processing ...")
                    plt.draw()
                    plt.pause(1/FRAMERATE) # framerate

                    # Slow edge colouring function
                    new_e_cols, new_num_cycles, new_e_freqs, endpoints = generate_edge_data(new_n, path_mode=path_mode)

                    # Get new graph data
                    new_faces = list(get_faces(new_n))
                    outer_face = max(new_faces, key=lambda x: perimeter(new_pos, x))
                    new_faces.remove(outer_face)
                    new_faces.append(outer_face)

                    # Update graph history
                    history = history[:history_index+1]
                    history.append((
                        deepcopy(new_n),
                        deepcopy(new_pos),
                        deepcopy(new_faces),
                        deepcopy(new_e_cols),
                        new_num_cycles,
                        deepcopy(endpoints),
                        deepcopy(new_e_freqs)
                    ))
                    history_index+=1
                    n, pos, faces, e_cols, num_cycles, endpoints, e_freqs = history[history_index]

                    # Export graph
                    if output is not None:
                        write_graph(pl_object, n)
                        output_index += 1

        EVENT_STATES['CLICK_EVENT'] = False

        # Change button colours if available
        BUTTONS['back'].is_valid = history_index != 0
        BUTTONS['forward'].is_valid = history_index != len(history)-1

        for b in BUTTONS.values():
            b.update_button()

        if EVENT_STATES['HOME_BUTTON']:
            EVENT_STATES['HOME_BUTTON'] = False
            history_index = 0
            n, pos, faces, e_cols, num_cycles, endpoints, e_freqs = history[history_index]

        # Handle redraw button event
        if BUTTONS['redraw'].state:
            BUTTONS['redraw'].set_state(False)
            g = nx.PlanarEmbedding()
            g.set_data(n)
            initial_pos = {k: np.array(v) for k, v in nx.planar_layout(g).items()}
            pos = scale_unit(initial_pos)

        # Handle back button event
        if BUTTONS['back'].state:
            BUTTONS['back'].set_state(False)
            if history_index > 0:
                history_index -= 1
                n, pos, faces, e_cols, num_cycles, endpoints, e_freqs = history[history_index]

        # Handle forward button event
        if BUTTONS['forward'].state:
            BUTTONS['forward'].set_state(False)
            if history_index < len(history) - 1:
                history_index += 1
                n, pos, faces, e_cols, num_cycles, endpoints, e_freqs = history[history_index]

        # Redraw edge and vertex positions
        v_obj.remove()
        for obj in t_objs.values():
            obj.remove()
        for obj in e_objs.values():
            obj.remove()
        for obj in e_text_objs.values():
            obj.remove()
        v_obj, t_objs, e_objs, e_text_objs = draw_from_coords(
            n,
            pos,
            colours=e_cols,
            e_freqs=e_freqs,
            num_cycles=num_cycles,
            do_endpoints=endpoints
        )

        #Reorder faces so the outer face is the last one
        outer_face = max(faces, key=lambda x: perimeter(pos, x))
        faces.remove(outer_face)
        faces.append(outer_face)

        # Update output flavour text if needed
        if output:
            save_message = f"Most recent graph is saved to {output.name}, index {output_index}"
        if path_mode == "L":
            path_name = f"Longest Path{plural(num_cycles)} with length {get_longest_path_length()}"
        elif path_mode == "H":
            path_name = f"Hamiltonian Cycle{plural(num_cycles)}"
        else:
            path_name = ""
        save_text.set_text(
            '' if path_mode == 'N' else f"Current graph has {num_cycles} {path_name}. {save_message}"
        )

        # Update back and forward button quantities
        BUTTONS['back'].set_text(f"Back ({history_index})")
        BUTTONS['forward'].set_text(f"Forward ({len(history) - 1 - history_index})")

        # Redraw all at a fixed framerate
        plt.draw()
        plt.pause(0.05)
    # Close active io workers
    plt.close()
    if output:
        pl_worker.__exit__()

def initiate(pl_name, graph_index=0, output=None, path_mode='H'):
    """Initiate the interactive loop"""
    graph_obj = read_graphs(pl_name)
    for _ in range(graph_index):
        try:
            next(graph_obj) # open correct graph
        except StopIteration:
            print(f"Error: file {pl_name} does not contain {'enough' if graph_index else 'any'} graphs.")
    graph_adj = next(graph_obj)

    while True:
        # Loop until plot closed
        print(graph_adj)
        draw_interactive(deepcopy(graph_adj), output=output, path_mode=path_mode)
        if not EVENT_STATES['HOME_BUTTON']:
            break
        EVENT_STATES['HOME_BUTTON'] = False
        print("Restarting")

def main(arguments):
    """Main for cli"""
    index = 0
    output = None
    path_mode = 'H'
    try:
        pl_name = Path(arguments[1])
        for i in range(1, len(arguments)):
            if arguments[i].lower() in ["--index", "-i"]:
                index = int(arguments[i+1])
            if arguments[i].lower() in ["--output", "-o"]:
                output = Path(arguments[i+1])
            if arguments[i].lower() in ["--longpaths", "-l"]:
                path_mode = 'L'
            if arguments[i].lower() in ["--nothing", "-n"]:
                path_mode = 'N'
    except IndexError:
        print(f"Usage: python3 {arguments[0]} <pl_name> [--index x] [--output filename] [--longpaths | --nothing]")
    else:
        initiate(pl_name, graph_index=index, output=output, path_mode=path_mode)

if __name__ == "__main__":
    main(argv)
