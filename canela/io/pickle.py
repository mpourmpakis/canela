import bz2
import pickle


def read_pickle_bz2(path):
    data = bz2.BZ2File(path, 'rb')
    data = pickle.load(data)
    return data


def read_pickle(path):
    with open(path, 'rb') as handle:
        data = pickle.load(handle)
    return data


def write_pickle(data, path):
    """
    Write arbitrary input data structure to disk as pickle.
    Enforces .pkl file extension
    """
    with open(path + '.pkl', 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)


def write_pickle_bz2(data, path):
    """
    Write arbitrary input data structure to disk as bz2 compressed pickle.
    Enforces .pkl.bz2 file extension
    """
    with bz2.BZ2File(path + '.pkl.bz2', 'wb') as handle:
        pickle.dump(data, handle)
