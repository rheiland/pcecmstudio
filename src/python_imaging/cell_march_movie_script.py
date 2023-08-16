import sys
import matplotlib.pyplot as plt
sys.path.append(r'../python_imaging')

from image_processing_for_physicell import *

options_for_figure2a = {}

options_for_figure2a = {"output_plot" : True,
                       "show_plot" : False,
                       "produce_for_panel" : False,
                        "plot_ECM_anisotropy" : True,
                        "plot_ECM_orientation" : True,
                        "retrieve_ECM_data": True,
                        "load_full_physicell_data" : True,
                        "plot_cells_from_SVG" : False,
                        "load_SVG_data": False,
                        'plot_cells_from_physicell_data': True,
                        "produce_for_movie" : True,
                        "contour_options" : {'lowest_contour': 0.90, 
                                           'upper_contour': 0.92,
                                           'number_of_levels': 25,
                                           'color_map_name': 'Reds',
                                           'color_bar': True
                                           },
                        "quiver_options" : {"scale_quiver": False,
                                          "mask_quiver": False}
                       }

movie_options_for_figure_2a = {}

movie_options_for_figure_2a = {'INCLUDE_ALL_SVGs': True,
                            'INCLUDE_FULL_HISTORY': False
                            }

#### Right now, if you don't have None or the full contour and quiver options, it will break in the plotting.

mf = PhysiCellPlotter()

mf.produce_movie(save_name='figure_2a_march', movie_options=movie_options_for_figure_2a, image_options=options_for_figure2a)

# generic_plotter (start, intervnal, finish, save_filename, data_path, save_path, options)
#
#     All based on options/logic- function
#     load_cell_positiondata
#     load_uE_data_chemical
#     load_uE_data_ECM
#
#     process data into plots - functions
#     - cell tracks (might be loaded by just be plotted???)
#     - cell positions
#     - ECM layer
#     - chemical layer
#
#     complete plot presentaiont and save (maybe functions)
#     - title
#     - axes