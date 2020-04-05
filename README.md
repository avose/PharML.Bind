# pharml
PharML is a framework for predicting interaction between protein structures and ligands. It utilizes a novel Molecular-Highway Graph Neural Network (MH-GNN) architecture based on state-of-the-art techniques in deep learning. This repository contains visualization, preprocessing, training, and inference code written in Python and C. In addition, it provides an ensemble of pre-trained models which can be used for generating predictions of compound binding relative to a given target.

Setup
==============================

1) Edit and execute the create_env_ubuntu.sh script considering your system configuration

    emacs tools/create_env_ubuntu.sh
    cd tools ; bash create_env_ubuntu.sh

3) Preprocess the Dataset(s)

   -> cd datasets/covid19 ; make

   -> After preprocessing completes, you will have a directory containing
        -> data/lig: Ligand graph files
        -> data/nhg: Protein neighborhood graph files
        -> data/pdb: Raw pdb files used to generate ligands / NHGs
        -> data/map: Map file used for inference, ligand-to-target tests

4) Test Inference Across example map file

    -> cd pharml.bind/ ; bash inference.sh
