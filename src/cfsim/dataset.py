import pandas as pd 
import h5py

class DataSet(object):

    @staticmethod
    def from_file(file_name):
        with h5py.File(file_name, "r") as fh:
            baf = fh["baf"][()]

            reads = fh["reads"][()]

            normal_reads = fh["normal"][()]

            bins = [x.decode() for x in fh["bins"]]

            cells = [x.decode() for x in fh["cells"]]

            gc = fh["gc"][()]

            mappability = fh["map"][()]

            regions = fh["region"][()]

            valid = fh["valid"][()]

        return DataSet(bins, cells, baf, reads, normal_reads, gc, mappability, regions, valid)

    def __init__(self, bins, cells, baf, reads, normal_reads, gc, mappability, regions, valid):
        self.bins = bins

        self.cells = cells

        self.baf = baf

        self.reads = reads

        self.normal_reads = normal_reads

        self.gc = gc

        self.map = mappability

        self.regions = regions

        self.valid = valid

    @property
    def bin_df(self):
        df = pd.DataFrame([x.split(":") for x in self.bins], columns=["chrom", "beg", "end"])

        df["chrom"] = df["chrom"].astype(str)

        df["beg"] = df["beg"].astype(int)

        df["end"] = df["end"].astype(int)

        return df

    @property
    def genome_size(self):
        df = self.bin_df

        df["len"] = df["end"] - df["beg"]

        return df["len"].sum()

    @property
    def num_bins(self):
        return self.baf.shape[1]

    @property
    def num_blocks(self):
        return self.baf.shape[2]

    @property
    def num_cells(self):
        return self.baf.shape[0]

    def filter_bins(self, bins):
        idxs = [self.bins.index(x) for x in bins]

        self.bins = [self.bins[i] for i in idxs]

        self.baf = self.baf[:, idxs]

        self.reads = self.reads[:, idxs]

        self.normal_reads = self.normal_reads[idxs]

        self.gc = self.gc[idxs]

        self.map = self.map[idxs]

        self.regions = self.regions[idxs]

        self.valid = self.valid[idxs]

    def filter_cells(self, cells):
        idxs = [self.cells.index(x) for x in cells]

        self.cells = [self.cells[i] for i in idxs]

        self.baf = self.baf[idxs]

        self.reads = self.reads[idxs]

    def remove_centromeres(self):
        bins = [self.bins[i] for i, r in enumerate(self.regions) if r != -1]

        self.filter_bins(bins)

