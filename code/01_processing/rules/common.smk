import os


def define_final_output(config, subset, out_path, suffix):
    config_subset = config['dataset'][subset]
    final_output = expand(
        os.path.join(out_path, config_subset['name_pattern'] + suffix),
        samplename = config_subset['samplenames'],
        rep = config_subset['reps'],
        exp = config_subset['exps'],
        tag = config_subset['tags'],
        num = [1, 2]
    )
    return final_output

def get_adapters(config, subset):
    config_subset = config['dataset'][subset]
    adapters = config['adapters'][config_subset['adapters']]
    return adapters

def get_refs(config, subset):
    config_subset = config['dataset'][subset]
    if isinstance(config_subset['ref'], dict):
        ref = dict([[x, config['refs'][y]] for x, y in config_subset['ref'].items()])
    else:
        ref = config['refs'][config_subset['ref']]
    return ref

def get_adapter_spikeins(config, subset):
    config_subset = config['dataset'][subset]
    if isinstance(config_subset['ref'], dict) and 'spikeins' in config_subset['ref'].keys():
        adapter_spikeins = config['adapters'][config_subset['ref']['spikeins']]
    else :
        adapter_spikeins = None
    return adapter_spikeins

def define_bowtie2_input(adapter_spikeins):
    if adapter_spikeins is None:
        return([os.path.join(config['dir']['results']['fq'], '{sample}_3end_1.fq.gz'), os.path.join(config['dir']['results']['fq'], '{sample}_3end_2.fq.gz')])
    else:
        return([os.path.join(config['dir']['results']['fq'], '{sample}_{tag}_1.fq.gz'), os.path.join(config['dir']['results']['fq'], '{sample}_{tag}_2.fq.gz')])
