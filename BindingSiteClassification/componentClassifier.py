import csv
import numpy as np
import plotly.figure_factory as ff
from plotly.offline import plot


def read_csv(file, heading=True):
    if heading:
        with open(file, 'r') as csv:
            lines = tuple(tuple(line.strip().split(',')) for line in csv.readlines()[1:])
        return lines
    else:
        with open(file, 'r') as csv:
            lines = tuple(tuple(line.strip().split(',')) for line in csv.readlines())
        return lines


csv_file = 'Z:\\project\\Metadata\\CSV\\best_weights_PIP_Entropy.csv'
csv_lines = read_csv(csv_file, heading=True)
AA_names = tuple(line[0] for line in csv_lines)
SM = tuple(tuple(map(float, line[1:4])) for line in csv_lines)
SM = np.array(SM)
fig = ff.create_dendrogram(X=SM, labels=list(AA_names), `color_threshold=0.5)
fig['layout'].update({'width': 1240})
fig['layout'].update({'height': 500})
ranking = tuple(fig['layout']['xaxis']['ticktext'])
AA_dict = {}
for line in csv_lines:
    AA_dict[line[0]] = line[1:]
rows = []
for aa in ranking:
    rows.append(AA_dict[aa])
rows.reverse()
"""with open('Z:\\project\\Metadata\\AA_clustering_vDW.csv', 'w+', newline='') as w:
    w = csv.writer(w)
    w.writerows(rows)"""
plot(fig)