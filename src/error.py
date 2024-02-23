"""Custom error classes."""


class InvalidSMILESError(Exception):
    """
    Exception raised when an invalid SMILES is encountered
    """


class InvalidSMARTSError(Exception):
    """
    Exception raised when Molecule can not be converted to SMARTS
    """