from itertools import permutations
from typing import Union
import numpy as np
import matplotlib.pyplot as plt
import ase
import canela.lpnc


def visualize_bond_distance_analysis(atoms: Union[str, ase.Atoms],
                                     save_path: str = None,
                                     show: bool = True,
                                     title: str = None,
                                     scale=canela.lpnc.SCALE) -> None:
    """visualizes a bond analysis based on atom ids
    - see canela.lpnc.get_atom_ids function for more details on
      id analysis.
    """
    # instantiate LPNC object
    mync = canela.lpnc.LPNC(atoms, scale=scale)

    # get the bond analysis dict
    bond_analysis_dict = canela.lpnc.get_bond_distance_analysis(mync.atoms,
                                                                scale=scale)

    # VISUALIZE
    fig, ax = plt.subplots(figsize=(16, 8))

    colors = {'MM': 'dodgerblue',
              'MS': 'lightgreen',
              'RS': 'orange',
              'MR': 'lightgray'}

    hatch = {'CC': '',
             'CS': '/',
             'SS': 'x'
             }

    names = {'C': 'Core',
             'S': 'Shell'
             }

    used_color = set()
    used_hatch = set()

    all_ids = set(mync.ids) - {'S_R_xx_xx'}
    x = 1
    xs = []
    labels = []
    used = set()
    for a, b in permutations(sorted(all_ids), 2):
        sortab = ''.join(sorted([a, b]))
        if sortab in used:
            continue

        if a in bond_analysis_dict and b in bond_analysis_dict[a]:
            y = np.array(bond_analysis_dict[a][b])

            # use M for metal if Core atom,
            # else use element symbol from id string
            ckey_list = sorted('M' if atom_id[0] == 'C'
                               else atom_id.split('_')[1]
                               for atom_id in [a, b])
            # join the symbols to create the color key,
            # which is also the label used in the legend
            ckey = ''.join(ckey_list)

            hkey = ''.join(sorted([a[0], b[0]]))
            used_color.add(ckey)
            used_hatch.add(hkey)

            ax.bar(x, y.mean(), yerr=y.std(), align='center',
                   color=colors[ckey], edgecolor='k', hatch=hatch[hkey])
            ax.scatter([x] * len(y), y, edgecolor='k', s=100, alpha=0.5,
                       color=colors[ckey], zorder=100)

            xs.append(x)
            asimp = a.strip('x').strip('_').strip('x').strip('_')
            bsimp = b.strip('x').strip('_').strip('x').strip('_')

            # create label to define bond type
            label = canela.lpnc.get_atom_id_latex_name(a)
            label += ' : '
            label += canela.lpnc.get_atom_id_latex_name(b)

            labels.append(label)
            x += 1
            used.add(sortab)
    ax.set_xticks(xs)
    ax.set_ylim(2, 3.6)
    ax.set_xticklabels(labels, rotation=60, fontsize=14)
    ax.set_ylabel('Bond Length $(\\AA)$', fontsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # add title if given
    if title:
        ax.set_title(title)

    for c in sorted(used_color):
        ax.bar(1, 0, color=colors[c], label=c, edgecolor='k')

    for h in sorted(used_hatch):
        ax.bar(1, 0, color='white', hatch=hatch[h] * 2,
               label='-'.join(names[s] for s in h), edgecolor='k')

    ax.legend(fontsize=14)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)
        print(f'Saved bond analysis to: {save_path}')
    if show:
        plt.show()
