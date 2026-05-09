import matplotlib.pyplot as plt
import matplotlib as mpl

def set_npj_style(column_type='single'):
    """
    Configures Matplotlib to meet Nature Portfolio (npj) standards.
    
    Args:
        column_type (str): 'single' (89mm) or 'double' (183mm).
    
    Returns:
        tuple: (width_in_inches, height_in_inches) to be passed to plt.figure()
    """
    
    # 1. Geometry (Nature Portfolio Standards)
    # 89 mm = 3.50 inches | 183 mm = 7.20 inches
    if column_type == 'double':
        fig_width_mm = 183
    else:
        fig_width_mm = 89
        
    fig_width_in = fig_width_mm / 25.4
    # Golden ratio is a good default for height, but adjustable
    fig_height_in = fig_width_in / 1.618 

    # 2. Typography
    # Nature requires sans-serif fonts (Arial or Helvetica preferred)
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans', 'Bitstream Vera Sans']
    mpl.rcParams['mathtext.default'] = 'regular' # Use same font for math as text
    mpl.rcParams['font.size'] = 8           # Base font size (increased to 8)
    mpl.rcParams['axes.labelsize'] = 9      # Increased to 9
    mpl.rcParams['axes.titlesize'] = 8      # Increased to 8
    mpl.rcParams['xtick.labelsize'] = 7     # Increased to 7
    mpl.rcParams['ytick.labelsize'] = 7     # Increased to 7
    mpl.rcParams['legend.fontsize'] = 6     # Decreased to 6 per user request
    mpl.rcParams['figure.titlesize'] = 9    # Increased to 9

    # 3. Lines and Markers
    mpl.rcParams['axes.linewidth'] = 0.5    # Edge line width
    mpl.rcParams['lines.linewidth'] = 1.0   # Data line width
    mpl.rcParams['lines.markersize'] = 3
    mpl.rcParams['xtick.major.width'] = 0.5
    mpl.rcParams['ytick.major.width'] = 0.5
    mpl.rcParams['xtick.minor.width'] = 0.4
    mpl.rcParams['ytick.minor.width'] = 0.4
    
    # 4. Layout & Aesthetics
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    
    # Systematically increase subplot spacing (constrained_layout padding)
    mpl.rcParams['figure.constrained_layout.use'] = True
    mpl.rcParams['figure.constrained_layout.w_pad'] = 0.1 # Increased horizontal spacing
    mpl.rcParams['figure.constrained_layout.h_pad'] = 0.1 # Increased vertical spacing
    
    mpl.rcParams['figure.dpi'] = 300         # High resolution for review
    mpl.rcParams['savefig.dpi'] = 600        # Publication quality
    mpl.rcParams['savefig.transparent'] = True
    mpl.rcParams['savefig.bbox'] = 'tight'
    mpl.rcParams['savefig.pad_inches'] = 0.05
    
    # 5. Legend
    mpl.rcParams['legend.frameon'] = False   # Cleaner look, preferred by Nature
    mpl.rcParams['legend.loc'] = 'best'
    
    return (fig_width_in, fig_height_in)

def add_panel_label(ax, label, x=-0.18, y=1.0):
    """
    Adds a panel label (a, b, c) in a consistent style and position.
    
    Args:
        ax (matplotlib.axes.Axes): The axes to label.
        label (str): The label text (e.g., 'a', 'b').
        x (float): X position in standard axes coordinates. Default -0.18.
        y (float): Y position in standard axes coordinates. Default 1.0.
    """
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='center')

def get_cpk_colors():
    """Returns a dictionary of standard CPK atom colors for computational materials."""
    return {
        'H': '#FFFFFF', 'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D',
        'F': '#90E050', 'Si': '#F0C8A0', 'P': '#FF8000', 'S': '#FFFF30',
        'Cl': '#1FF01F', 'Br': '#A62929', 'I': '#940094', 'Li': '#CC80FF', 
        'Na': '#AB5CF2', 'K': '#8F40D4', 'Fe': '#E06633', 'Cu': '#C88033'
    }