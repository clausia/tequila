"""
Export QCircuits as qasm code

OPENQASM version 2.0 specification from:
A. W. Cross, L. S. Bishop, J. A. Smolin, and J. M. Gambetta, e-print arXiv:1707.03429v2 [quant-ph] (2017).
https://arxiv.org/pdf/1707.03429v2.pdf
"""
from tequila import TequilaException
from tequila.circuit import QCircuit
from tequila.circuit.compiler import Compiler
from tequila.circuit.gates import *
import re


def export_open_qasm(circuit: QCircuit, variables=None, version: str = "2.0", filename: str = None, zx_calculus: bool = False) -> str:
    """
    Allow export to different versions of OpenQASM

    Args:
        circuit: to be exported to OpenQASM
        variables: optional dictionary with values for variables
        version: of the OpenQASM specification, optional
        filename: optional file name to save the generated OpenQASM code
        zx_calculus: indicate if y-gates must be transformed to xz equivalents

    Returns:
        str: OpenQASM string
    """

    if version == "2.0":
        result = convert_to_open_qasm_2(circuit=circuit, variables=variables, zx_calculus=zx_calculus)
    else:
        return "Unsupported OpenQASM version : " + version
    # TODO: export to version 3

    if filename is not None:
        with open(filename, "w") as file:
            file.write(result)
            file.close()

    return result


def import_open_qasm(qasm_code: str, variables=None, version: str = "2.0", rigorous: bool = True) -> QCircuit:
    """
    Allow import from different versions of OpenQASM

    Args:
        qasm_code: string with the OpenQASM code
        variables: optional dictionary with values for variables
        version: of the OpenQASM specification, optional
        rigorous: indicates whether the QASM code should be read rigorously

    Returns:
        QCircuit: equivalent to the OpenQASM code received
    """

    # print("Import from Open QASM", qasm_code)

    lines = qasm_code.splitlines
    oq_code = []
    # ignore comments
    for line in lines:
        if line.find("//") != -1:
            clean_line = line[0:line.find("//")].strip()
        else:
            clean_line = line.strip()
        if clean_line:
            oq_code.append(clean_line)

    if oq_code[0].startswith("OPENQASM"):
        oq_code.pop(0)
    elif rigorous:
        raise TequilaException("File must start with the 'OPENQASM' directive")

    if oq_code[0].startswith('include "qelib1.inc";'):
        oq_code.pop(0)
    elif rigorous:
        raise TequilaException("File must import standard library")

    code = "\n".join(oq_code)
    # separate the custom command definitions from the normal commands
    while True:
        i = code.find("gate ")
        if i == -1:
            break
        j = code.find("}", i)
        parse_custom_gate(code[i:j + 1])
        code = code[:i] + code[j + 1:]

    # parse regular commands
    commands = [s.strip() for s in code.split(";") if s.strip()]
    circuit = QCircuit()
    for c in commands:
        res = parse_command(c)
        if res is not None:
            circuit += res

    return circuit


def import_open_qasm_from_file(filename: str, variables=None, version: str = " 2.0", rigorous: bool = True) -> QCircuit:
    """
    Allow import from different versions of OpenQASM from a file

    Args:
        filename: string with the file name with the OpenQASM code
        variables: optional dictionary with values for variables
        version: of the OpenQASM specification, optional
        rigorous: indicates whether the QASM code should be read rigorously

    Returns:
        QCircuit: equivalent to the OpenQASM code received
    """

    with open(filename, "r") as file:
        qasm_code = file.read()
        file.close()

    return import_open_qasm(qasm_code, variables=variables, version=version, rigorous=rigorous)


def convert_to_open_qasm_2(circuit: QCircuit, variables=None, zx_calculus: bool = False) -> str:
    """
    Allow export to OpenQASM version 2.0

    Args:
        circuit: to be exported to OpenQASM
        variables: optional dictionary with values for variables
        zx_calculus: indicate if y-gates must be transformed to xz equivalents

    Returns:
        str: OpenQASM string
    """

    if variables is None and not (len(circuit.extract_variables()) == 0):
        raise TequilaException(
            "You called export_open_qasm for a parametrized type but forgot to pass down the variables: {}".format(
                circuit.extract_variables()))

    compiler = Compiler(multitarget=True,
                        multicontrol=False,
                        trotterized=True,
                        generalized_rotation=True,
                        exponential_pauli=True,
                        controlled_exponential_pauli=True,
                        hadamard_power=True,
                        controlled_power=True,
                        power=True,
                        toffoli=True,
                        controlled_phase=True,
                        phase=True,
                        phase_to_z=True,
                        controlled_rotation=True,
                        swap=True,
                        cc_max=True,
                        gradient_mode=False,
                        ry_gate=zx_calculus,
                        y_gate=zx_calculus)

    compiled = compiler(circuit, variables=None)

    result = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n"

    qubits_names = dict()
    for q in compiled.qubits:
        name = "q[" + str(q) + "]"
        qubits_names[q] = name

    result += "qreg q[" + str(compiled.n_qubits) + "];\n"
    result += "creg c[" + str(compiled.n_qubits) + "];\n"

    for g in compiled.gates:

        control_str = ''
        if g.is_controlled():

            if len(g.control) > 2:
                raise TequilaException(
                    "Multi-controls beyond 2 not yet supported for OpenQASM 2.0. Gate was:\n{}".format(g))

            controls = list(map(lambda c: qubits_names[c], g.control))
            control_str = ','.join(controls) + ','

        gate_name = name_and_params(g, variables)
        for t in g.target:
            result += gate_name
            result += control_str
            result += qubits_names[t]
            result += ";\n"

    return result


def name_and_params(g, variables):
    """
    Determines the quantum gate name and its parameters if applicable

    Args:
        g: gate to get its name
        variables: dictionary with values for variables

    Returns:
        str: name (and parameter) to the gate specified
    """

    res = ""

    for c in range(len(g.control)):
        res += "c"

    res += g.name.lower()

    if hasattr(g, "parameter") and g.parameter is not None:
        res += "(" + str(g.parameter(variables)) + ")"

    res += " "

    return res


def parse_custom_gate(custom: str) -> None:
    """
    Parse custom gates code

    Args:
        custom: code with custom gates
    """
    len(custom)


def parse_command(command: str) -> QCircuit:
    """
    Parse qasm code

    Args:
        command: open qasm code to be parsed
    """

    name, rest = command.split(" ", 1)
    if name in ("barrier", "creg", "measure", "id"):
        return None
    if name in ("opaque", "if"):
        raise TequilaException("Unsupported operation {}".format(command))
    args = [s.strip() for s in rest.split(",") if s.strip()]

    # if name == "qreg":
    #     regname, sizep = args[0].split("[", 1)
    #     size = int(sizep[:-1])
    #     registers[regname] = (qubit_count, size)
    #     qubit_count += size
    #     return

    if name == "x":
        return X(target=int(re.search(r"\[([0-9]+)\]", args[0]).group(1)))
    if name == "y":
        return Y(target=int(re.search(r"\[([0-9]+)\]", args[0]).group(1)))
    if name == "z":
        return Z(target=int(re.search(r"\[([0-9]+)\]", args[0]).group(1)))
