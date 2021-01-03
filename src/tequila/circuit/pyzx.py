"""
Add to tequila the ability to make ZX-Calculus

Using the pyzx library: https://github.com/Quantomatic/pyzx
"""

HAS_PYZX = True
try:
    import pyzx
    HAS_PYZX = True
except ImportError:
    HAS_PYZX = False


from tequila import export_open_qasm, import_open_qasm
from tequila.circuit import QCircuit


def convert_to_pyzx(circuit: QCircuit, variables=None) -> pyzx.circuit.Circuit:
    """
    Allow convert from Tequila circuit to pyzx circuit

    Args:
        circuit: in Tequila format to be exported to pyzx
        variables: optional dictionary with values for variables

    Returns:
        pyzx.circuit.Circuit: pyzx circuit
    """
    if HAS_PYZX:
        return pyzx.circuit.Circuit.from_qasm(export_open_qasm(circuit=circuit, variables=variables, version="2.0"))
    else:
        raise TequilaException("Pyzx package not installed.")


def convert_from_pyzx(circuit: pyzx.circuit.Circuit) -> QCircuit:
    """
    Allow convert from pyzx circuit to Tequila circuit

    Args:
        circuit: in pyzx format to be exported to Tequila circuit

    Returns:
        QCircuit: Tequila circuit
    """
    if HAS_PYZX:
        return import_open_qasm(circuit.to_qasm(), version="2.0")
    else:
        raise TequilaException("Pyzx package not installed.")