"""The cnstools package can be run in three ways. It can be run as an executable directly from the commandline, 
run as a python program with ``$python cnstools`` or it can be imported as a python module :class:`cnstools` and used with other 
python scripts. Currently this documentation only covers the imported module, more to come later."""
import align, score, identify, analyze
__all__ = [
    "align",
    "score",
    "identify",
    "analyze"
]