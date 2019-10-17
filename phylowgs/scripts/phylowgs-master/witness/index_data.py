from __future__ import print_function
import os
import glob
import json
import os
from collections import defaultdict

def main():
  datasets = defaultdict(list)

  base_dir = 'data'
  for run_name in os.listdir(base_dir):
    for summary_path in glob.glob(os.path.join(base_dir, run_name, '*.summ.json')):
      dataset_name = summary_path.split('/')[-1].split('.')[0]
      run = {
        'summary_path': summary_path,
        'muts_path': os.path.join(base_dir, run_name, dataset_name + '.muts.json'),
        'name': dataset_name,
      }
      print(run)
      clusters_path = os.path.join(base_dir, run_name, dataset_name + '.clusters.json')
      if os.path.exists(clusters_path):
        run['clusters_path'] = clusters_path
      datasets[run_name].append(run)

  out_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'index.json')
  with open(out_path, 'w') as outf:
    print(json.dumps(datasets), file=outf)

main()
