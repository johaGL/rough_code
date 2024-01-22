# copied from the official pyCirclize github
from pycirclize import Circos
import pandas as pd
import numpy as np

# Create matrix dataframe (3 x 6)
row_names = ["T1", "T2", "T3", "T4", "T5", "T6", "7", "8", "9", "10", "11", "12", "13"]
col_names = ["DE:-log10 padj" ]

# attention: this does not accept values of zero, so -log10 padj when zero--> 0.1
matrix_data = 	[[  22],
	 [   20],
	 [  15],
	 [  10],
	 [   0.1],
	 [   0.1],
	 [  0.1],
	 [  0.1],
	 [   0.1],
	 [  20],
	 [  22],
	 [  24],
	 [ 24]]



matrix_df = pd.DataFrame(matrix_data, index=row_names, columns=col_names)

colors_dict = { "T1" : "darkorchid",  "T2": "salmon", "T3": "lightcoral", "T4": "bisque", "T5": "white", 
    "T6": "white", "7": "white", "8": "white", "9": "white", "10": "white", "11": "lightblue",
    "12": "royalblue", "13": "blue"

}
dict(A="royalblue", B="orange", C="green", D="red"),
# Initialize Circos from matrix for plotting Chord Diagram
circos = Circos.initialize_from_matrix(
    matrix_df,
    space=3,
    #cmap="tab10",
    cmap=colors_dict,
    r_lim=(93, 100),
    label_kws=dict(size=12),
    link_kws=dict(ec="black", lw=0.5, direction=0),
)


circos.savefig("chord3.svg")
