# Chrysanthemic Acid Biosynthesis Model

This project includes all the scripts and files required for the model simulations in the paper: **"De novo biosynthesis of chrysanthemic acid in *Escherichia coli* with combinational metabolic modeling and engineering"**.

## Project Structure

### `Model` Folder
Contains the original model (pre-reconstruction) and the reconstructed new model.

### `Database` Folder
Contains essential data for running the scripts, including:
- **Genome CDS sequences** (used for manual curation of the model),
- **BLASTp results** (used for gap filling in the model),
- **Model database files** (used to correct metabolite and reaction information in the model).

### `Memote Report` Folder
Stores the memote reports of the model, which provide an evaluation of the model's quality and integrity.

### `Script` Folder
Contains the functions and scripts used throughout the project for model simulations and data processing tasks.
