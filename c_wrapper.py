import ctypes
from numpy.ctypeslib import ndpointer
c_lib = ctypes.cdll.LoadLibrary('./reresquiggle.so')

new_KmerModel = c_lib.new_KmerModel
new_KmerModel.restype = ctypes.c_void_p
new_KmerModel.argtypes = [ctypes.c_int,
                         ctypes.c_int,
                         ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                         ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')]

delete_KmerModel = c_lib.delete_KmerModel
delete_KmerModel.restype = None
delete_KmerModel.argtypes = [ctypes.c_void_p]

new_ResquiggledRead = c_lib.new_ResquiggledRead
new_ResquiggledRead.restype = ctypes.c_void_p
new_ResquiggledRead.argtypes = [ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                               ctypes.c_int,
                               ndpointer(ctypes.c_int, flags='C_CONTIGUOUS'),
                               ndpointer(ctypes.c_int, flags='C_CONTIGUOUS'),
                               ctypes.c_int,
                               ctypes.c_int]

delete_ResquiggledRead = c_lib.delete_ResquiggledRead
delete_ResquiggledRead.restype = None
delete_ResquiggledRead.argtypes = [ctypes.c_void_p]

compute_probabilities = c_lib.compute_probabilities
compute_probabilities.restype = None
compute_probabilities.argtypes = [ndpointer(ctypes.c_int, flags='C_CONTIGUOUS'),
                                  ctypes.c_int,
                                  ctypes.c_void_p,
                                  ctypes.c_void_p,
                                  ndpointer(ctypes.c_int, flags='C_CONTIGUOUS'),
                                  ctypes.c_void_p,
                                  ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_double,
                                  ctypes.c_bool]
