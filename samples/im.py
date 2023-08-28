import subprocess
import os
from graphs.py import graphs

k=100
R=256
if __name__ == '__main__':
    IM = f'{CURRENT_DIR}/../benchmark'
    out_file = f"{CURRENT_DIR}/../time.txt"
    for graph, url, w, GRAPH_DIR in graphs:
        graph_file = f'{GRAPH_DIR}/{graph}.bin'
        if not os.path.exists(graph_file):
            print(f'\nWarning: {graph} does not exists')
            continue
        print(f'\nRunning IM on {graph_file}')
        subprocess.call(f"echo '{graph}' >> {out_file}")
        cmd = f'{IM} {graph_file} {k} {R} -UIC {w} >> {out_file}'
        subprocess.call(cmd, shell=True)