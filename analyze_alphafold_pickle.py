#!/usr/bin/env python3
"""Analyze AlphaFold Pickle file and retrieve model characteristics.
"""

import pickle
import sys
import os
import argparse
import matplotlib
import matplotlib.pyplot as plt


class ModelData:
    def __init__(self, pkl_file_name, multimer):
        "Parse necessary data from AlphaFold pickle file"
        self.file_name = os.path.splitext(pkl_file_name)[0]
        p = pickle.load(open(pkl_file_name, 'rb'))
        self.ptm = float(p['ptm'])
        if multimer:
            self.iptm = float(p['iptm'])
        self.plddt = p['plddt']
        self.global_plddt = sum(p['plddt'])/len(p['plddt'])
        self.predicted_aligned_error = p['predicted_aligned_error']
        self.max_pae = p['max_predicted_aligned_error']

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


def main():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('pkl_file')
    arguments_parser.add_argument(
        '--multimer', action='store_true', help='parse multimer results'
        )
    arguments_parser.add_argument(
        '--save-plots', action='store_true',
        help='save pLDDT and PAE plots to files'
        )
    args = arguments_parser.parse_args()

    model_data = ModelData(args.pkl_file, args.multimer)
    model_data.plot_PAE(args.save_plots)
    model_data.plot_pLDDT(args.save_plots)
    print(model_data)


if __name__ == '__main__':
    sys.exit(main())
