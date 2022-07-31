"""
Assemble and sort some COVID reads...
"""

import subprocess
from pathlib import Path
import os
import pickle as pickle
import sqlite3
import pandas as pd

from latch import small_task, workflow
from latch.types import LatchFile

from enum import Enum
from typing import Optional, List

class pfam_option(Enum):
    G_protein_family = 'G_protein_family'
    Beta_lactamase_family = 'Beta_lactamase_family'
    Other_family = 'Other_family'

class chainID_option(Enum):
    A = 'A'
    B = 'B'
    C = 'C'
    D = 'D'
    E = 'E'

@small_task
def sca_compute_task(
         data_fam_dir: Optional[LatchFile],
         data_pdb_dir: Optional[LatchFile],
         output_dir: Optional[str],
         family: pfam_option = pfam_option.G_protein_family,
         chainID: chainID_option = chainID_option.A,
         ) -> (LatchFile, LatchFile):


    # set output file
    if not output_dir:
        output_path = Path("").resolve()

    else:
        output_path = output_dir


    if str(family.value) == 'G_protein_family':
        _run_SCAprocess = [  "/root/pySCA_temp/bin/scaProcessMSA",
            "-a",
            "/root/pySCA_temp/data/PF00071_rd2.an",
            "-b",
            "/root/pySCA_temp/data",
            "-s",
            "5P21",
            "-c",
            "A",
            "-f",
            "'Homo sapiens'",
            "-t",
            "-n",
            "-d",
            "./output"
        ]

        _run_ANNOTprocess = None
        filename = "PF00071_rd2"
    
    elif str(family.value) == 'Beta_lactamase_family':
        _run_SCAprocess = [   "/root/pySCA_temp/bin/scaProcessMSA",
            "-a",
            "/root/pySCA_temp/data/PF13354_full.an",
            "-b",
            "/root/pySCA_temp/data",
            "-s",
            "1FQG",
            "-c",
            "A",
            "-f",
            "'Escherichia coli'",
            "-t",
            "-n",
            "-d",
            "./output"
        ]

        _run_ANNOTprocess = None
        filename = "PF13354_full"
        
    elif str(family.value) == 'Other_family':
        

        # latch file path
        load_data_path = data_fam_dir.remote_path
        # filename of the .an data
        load_filename = load_data_path.split('/')[-1]
        # filename
        filename = load_filename.replace('.an', '')
  
        # copy the sequences from latch onto /root/ 
        subprocess.run([
              "latch",
              'cp',
              load_data_path,
              f"/root/pySCA_temp/data/{load_filename}"
        ])

        # pdb structure of the protein of interest
        load_pdb_path = data_pdb_dir.remote_path
        # filename of the .pdb data
        load_pdb_filename = load_pdb_path.split('/')[-1]
        # pdb filename
        pdb_filename = load_pdb_filename.replace('.pdb', '')
         
         # copy the sequences from latch onto /root/ 
        subprocess.run([
              "latch",
              'cp',
              load_pdb_path,
              f"/root/pySCA_temp/data/{pdb_filename}"
        ])
   
        
        _run_SCAprocess = [  
            "/root/pySCA_temp/bin/scaProcessMSA",
            "-a",
            f"/root/pySCA_temp/data/{load_filename}",
            "-b",
            "/root/pySCA_temp/data",
            "-s",
            f"{pdb_filename}",
            "-c",
            chainID.value,
            "-t",
            "-n",
            "-d",
            "./output"
        ]

        _run_ANNOTprocess = None
        

    else:
        pass
    # move into the SCA repo:
    os.chdir("./pySCA_temp/bin")
    os.mkdir("./output")

    if _run_ANNOTprocess is None:
        pass
    else:
        subprocess.run(_run_ANNOTprocess)

    subprocess.run(_run_SCAprocess)

    _run_MOVEoutput = [
            "cp",
            "-r",
            "./output",
            "./bin/"
            ]
    subprocess.run(_run_MOVEoutput)

  
    output_SCAcore = "/root/pySCA_temp/bin/output/" + filename + ".db"
    _run_SCAcore = [
            "/root/pySCA_temp/bin/scaCore",
            "-i",
            output_SCAcore
    ]

    subprocess.run(_run_SCAcore)


    output_SCAsector = "/root/pySCA_temp/bin/output/" + filename + ".db"
    _run_SCAsector = [
            "/root/pySCA_temp/bin/scaSectorID",
            "-i",
            output_SCAsector
    ]
 
    subprocess.run(_run_SCAsector)



    return (LatchFile(f"/root/pySCA_temp/bin/output/{filename}.db", f'latch:///{output_path}{filename}.db'),
            LatchFile(f"/root/pySCA_temp/bin/output/{filename}_processed.fasta", f'latch:///{output_path}{filename}_processed.fasta')
    )

@workflow
def sca_workflow(
        data_fam_dir: Optional[LatchFile],
        data_pdb_dir: Optional[LatchFile],
        output_dir: Optional[str],
        family: pfam_option = pfam_option.G_protein_family,
        chainID: chainID_option = chainID_option.A,
        ) -> (LatchFile, LatchFile):
    """Statistical Coupling Analysis (SCA) is a statistical method for defining and extracting the connective pattern of coevolving amino acids in a protein family.

    # Statistical Coupling Analysis (SCA)
    This package provides an implementation of the inference pipeline of Statistical Coupling Analysis (SCA) by the Ranganathan lab. This is a statistical method that defines and extracts the connective pattern of coevolving amino acids in a protein family. This algorithm returns the statistical structure of the multiple sequence alignment (MSA), SCA conserviation and evolution, and identify sectors. For more details please check out the paper by the Ranganathan group: ["Evolution-Based Functional Decomposition of Proteins" (Olivier R. et al., 2016)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004817). The four pre-define protein families are great starting points: G-coupled receptor family, S1A serine protease family, Dihydrofolate reductase (DHFR) family, and Beta-Lactamase enzyme family. 


    ## Citations

    If you use the code or data in this package, please cite:

    ```

    @article{{rivoire2016evolution},
        title={Evolution-based functional decomposition of proteins},
        author={Rivoire, Olivier and Reynolds, Kimberly A and Ranganathan, Rama},
        journal={PLoS computational biology},
        volume={12},
        number={6},
        pages={e1004817},
        year={2016},
        publisher={Public Library of Science San Francisco, CA USA}
    }
  
    ```

    ## Liscense
    The original SCA source code is free software distributed under the BSD 3-clause license, please see the file LICENSE for details.

    Copyright (C) 2019 Olivier Rivoire, Rama Ranganathan, and Kimberly Reynolds.


    __metadata__:
        display_name: Statistical coupling analysis (SCA)
        author:
            name: Ranganathan Lab
            email:
            github: https://github.com/ranganathanlab/pySCA
        repository:
        license:
            id: BSD 3-clause liscense

    Args:

        data_fam_dir:
          (Optional) Data directory that contains the protein family of interest..

          __metadata__:
            display_name: Choose your protein family directory (optional).
        
        data_pdb_dir:
          (Optional) PDB data direcotry that contains the reference sequence's pdb structure

          __metadata__:
            display_name: Choose your pdb structure for the reference sequence

        output_dir:
          Output path containing the SCA results.

          __metadata__:
            display_name: Output directory containing the SCA results and figures.

        family:
          Protein family option.

          __metadata__:
            display_name: Pick between the four protein family examples or your own protein family of interest.
 
        chainID:
          chain ID in the PDB for the reference sequence.

          __metadata__:
            display_name: Choose your chain ID in the PDB for the reference sequence.

    """
    return sca_compute_task(
            data_fam_dir = data_fam_dir,
            data_pdb_dir = data_pdb_dir,
            output_dir = output_dir,
            family = family,
            chainID = chainID,
    )
