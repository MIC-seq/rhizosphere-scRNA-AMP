"""Project-local PhyloPhlAn launcher; preserve whitespace in path arguments."""
import time
from pathlib import Path
import phylophlan.phylophlan as engine

original_compose = engine.compose_command


def compose_checked(params, check=False, sub_mod=None, input_file=None,
                    database=None, output_path=None, output_file=None, nproc=1):
    values = dict(input_file=input_file, database=database,
                  output_path=output_path, output_file=output_file)
    replacements = {}
    for name, value in values.items():
        if value is not None and any(c.isspace() for c in str(value)):
            token = '__PHYLO_PATH_' + name.upper() + '__'
            replacements[token] = str(value)
            values[name] = token
    result = original_compose(params, check=check, sub_mod=sub_mod, nproc=nproc, **values)

    def restore(value):
        if isinstance(value, str):
            for token, path in replacements.items():
                value = value.replace(token, path)
        elif isinstance(value, list):
            value = [restore(item) for item in value]
        return value

    return {key: restore(value) for key, value in result.items()}


engine.compose_command = compose_checked
original_build = engine.build_phylogeny


def build_checked(configs, key, inputt, output_path, output_tree, nproc=1, verbose=False):
    prefix = Path(output_path) / output_tree
    tree = Path(str(prefix) + '.treefile')
    logfile = Path(str(prefix) + '.log')
    if tree.is_file() and logfile.is_file():
        log = logfile.read_text(errors='replace')
        if tree.stat().st_size and 'Total wall-clock time used:' in log:
            engine.info('Reusing successfully completed IQ-TREE output: ' + str(tree) + '\n')
            return
    return original_build(configs, key, inputt, output_path, output_tree,
                          nproc=nproc, verbose=verbose)


engine.build_phylogeny = build_checked

if __name__ == '__main__':
    start = time.time()
    engine.phylophlan_main()
    engine.info('\nChecked launcher completed in %.1f seconds\n' % (time.time() - start))
