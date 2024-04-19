#!/usr/bin/env python3
"""Analyze AlphaFold Pickle file and retrieve model characteristics.
"""

import pickle
import sys
import os
import argparse
import json
import numpy
import matplotlib
import matplotlib.pyplot as plt


class ModelData:
    def __init__(self, pkl_file_name, multimer):
        "Parse necessary data from AlphaFold pickle file"
        self.file_name = os.path.splitext(pkl_file_name)[0]
        p = pickle.load(open(pkl_file_name, 'rb'))
        self.multimer = multimer
        del p['distogram']
        del p['experimentally_resolved']
        del p['masked_msa']
        del p['structure_module']
        del p['aligned_confidence_probs']
        self.data = p

    @property
    def ptm(self):
        return float(self.data['ptm'])

    @property
    def iptm(self):
        if self.multimer:
            iptm = float(self.data['iptm'])
        else:
            iptm = None
        return iptm

    @property
    def plddt(self):
        return self.data['plddt']

    @property
    def global_plddt(self):
        return sum(self.data['plddt'])/len(self.data['plddt'])

    @property
    def predicted_aligned_error(self):
        return self.data['predicted_aligned_error']

    @property
    def max_pae(self):
        return float(self.data['max_predicted_aligned_error'])

    def __str__(self):
        s = []
        s.append('pLDDT %s' % self.global_plddt)
        s.append('pTM %s' % self.ptm)
        if hasattr(self, 'iptm'):
            s.append('ipTM %s' % self.iptm)
        return '\n'.join(s)

    def plot_PAE(self, save_to_file=False):
        plt.title('Predicted aligned error')
        plt.imshow(
            self.predicted_aligned_error,
            vmin=0, vmax=self.max_pae, cmap='bwr'
            )
        if save_to_file:
            plt.savefig(self.file_name+'_PAE.png')
        else:
            plt.show()
        plt.close()

    def plot_pLDDT(self, save_to_file=False):
        plt.title('Predicted LDDT')
        plt.ylim(0, 100)
        plt.plot(self.plddt)
        if save_to_file:
            plt.savefig(self.file_name+'_pLDDT.png')
        else:
            plt.show()
        plt.close()

    def to_json(self):
        "Convert AlphaFold pkl to json with necessary data only"
        return json.dumps(self.data, cls=NumPyEncoder)


class NumPyEncoder(json.JSONEncoder):
    def default(self, o):
        if type(o).__module__ == numpy.__name__:
            if isinstance(o, numpy.ndarray):
                return o.tolist()
            else:
                return o.item()
        else:
            return json.JSONEncoder.default(self, o)


def main():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('pkl_file')
    arguments_parser.add_argument(
        '--multimer', action='store_true', help='parse multimer results'
        )
    arguments_parser.add_argument(
        '--show-plots', action='store_true',
        help='show pLDDT and PAE plots'
        )
    arguments_parser.add_argument(
        '--save-plots', action='store_true',
        help='save pLDDT and PAE plots to files'
        )
    arguments_parser.add_argument(
        '--print-json', action='store_true',
        help='print model data in JSON format'
        )
    args = arguments_parser.parse_args()

    model_data = ModelData(args.pkl_file, args.multimer)
    if args.show_plots or args.save_plots:
        model_data.plot_PAE(args.save_plots)
        model_data.plot_pLDDT(args.save_plots)
    if args.print_json:
        print(model_data.to_json())
    else:
        print(model_data)


if __name__ == '__main__':
    sys.exit(main())
