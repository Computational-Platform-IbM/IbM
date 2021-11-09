import markdown
import io

def readmd():
    """Read simulation_log.md and return the last simulation number.

    Args:
        No args

    Returns:
        string: number of the last simulation
    """

    with open('simulation_log.md', 'r', encoding='utf-8') as f:
        ptext = io.StringIO(markdown.markdown(f.read()))
        lastline = ptext.readlines()[-1] # Last line
        matID = [int(s) for s in lastline.split() if s.isdigit()]
        f.close()
        
    return(matID[0])
    
    
def writemd(simID: int, simname:str, siminfo: str, simgoal: str, version: str, replicate: int):
    """Write the new simulation in simulation_log.md

    Args:
        simID (int): numberID of the new simulation
        simname (str): simulation name
        siminfo (str): simulation description
        simgoal (str): aim of the simulation
        version (str): version of IbM in which simulation is performed
        replicate (int): replicate namber
    """
    
    # newtext = | Simulation number | Simulation description | Goal | Finished (default: empty) | Version |
    newtext = '| {:0>4.0f} | <b>{}</b> (Run {:.0f}){} | {} | | {} |\n'.format(simID, simname, replicate, siminfo, simgoal, version)
    
    with open('simulation_log.md', 'a', encoding='utf-8') as f:
        f.write(newtext)
        f.close()