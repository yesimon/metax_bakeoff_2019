import copy
import os
import yaml
def merge_yaml_dicts(*dicts):
    """Merge yaml dicts.

    Only handle top-level yaml dicts. Simply uses dict.update and ignores None
    dicts.
    """
    out = copy.deepcopy(dicts[0])
    for d in dicts[1:]:
        if d:
            out.update(d)
            return out


def load_config(yaml_path):
    with open(yaml_path) as f:
        return yaml.load(f.read())

user_config_fn = os.path.expanduser('~/.config/snakemake/config.yaml')
if os.path.exists(user_config_fn):
    user_config = load_config(user_config_fn)
    local_config = config
    config = merge_yaml_dicts(user_config, local_config)


import multiprocessing
ALL_CORES = config.get('ALL_CORES', multiprocessing.cpu_count())
