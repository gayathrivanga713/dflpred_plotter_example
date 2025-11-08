What is DFLpred?

DFLpred is a computational prediction tool designed to identify disordered flexible linker (DFL) regions in protein sequences. These are regions that:

are intrinsically disordered (lack a stable 3D structure under physiological conditions),

function as linkers between structured domains or segments of a protein (rather than being tightly folded themselves), and

often provide flexibility, enable domain movement, or facilitate functional transitions.

Why it is useful?
Many proteins, especially multi‐domain ones, need linker regions that are flexible and disordered for proper function (e.g., enabling the relative motion of domains, or acting as flexible spacers). Using a tool like DFLpred helps identify likely linker sequences from sequence alone.

It complements other predictors of intrinsic disorder and functional regions (e.g., binding regions). In fact, DFLpred is one of the “disorder‐function” predictors (linker prediction) included in the DEPICTER2 toolkit. 

Helps in annotation of proteins (especially uncharacterised ones) by flagging segments likely to be flexible linkers, which can be helpful for modelling, mutagenesis design, or understanding domain organisation.

How it works?

Input: protein sequence(s) in FASTA format.

Feature generation: From the input sequence, various features are computed (profiles of sequence composition, predicted disorder, predicted secondary structure, perhaps conservation etc). (The exact details of all features used by DFLpred in the original paper.)

Prediction: Based on these features, a machine‐learning model (in the original DFLpred likely logistic/regression or similar) predicts, for each residue: a propensity score for being a disordered flexible linker residue and a binary classification (linker vs non-linker) (when thresholded).

Output: For each residue in the sequence you get a linker‐propensity value.

DFLpred Plotter:

This script plots **DFLpred** per-residue scores along a protein sequence, calls
low-flexibility (DFL) regions using a score threshold, and (optionally) compares them to
a simple **case-based annotation** convention (lowercase letters in the input sequence
are treated as DFL by annotation).

Requirements
- Python 3.8+
- matplotlib

bash
pip install matplotlib
