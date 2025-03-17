import numpy as np
import os

def itol_colorstrip(legend, annotation, sample_name, colors_dict, output_path):
    # Write itol file

    for value, df in legend.groupby(annotation):
        legend.loc[df.index, 'color'] = colors_dict[value]
    
    legend_labels = list(colors_dict.keys())
    colors = list(colors_dict.values())

    with open(output_path, "w") as itol_file:
        itol_file.write("DATASET_COLORSTRIP\n\n")
        itol_file.write("SEPARATOR TAB\n\n")
        itol_file.write("DATASET_LABEL\t{0}\nCOLOR\t{1}\n\n".format(annotation, list(colors_dict.values())[-1]))
        itol_file.write("LEGEND_TITLE\t{0}\nLEGEND_SHAPES\t{1}\nLEGEND_COLORS\t{2}\nLEGEND_LABELS\t{3}\n\n".format(annotation,
                                                                                                                   "\t".join(['1']*len(np.unique(legend[annotation]))),
                                                                                                                    "\t".join(colors),
                                                                                                                    "\t".join(str(x) for x in legend_labels)))
        itol_file.write("BORDER_WIDTH\t0.25\nBORDER_COLOR\t#CCCCCC\n\n")
        itol_file.write("DATA\n")
        for i,row in legend.iterrows():
            itol_file.write("{0}\t{1}\t{2}\n".format(row[sample_name], row['color'], row[annotation]))