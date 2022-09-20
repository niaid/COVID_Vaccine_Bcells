
Download all_subjects_cd19_positive_and_keys.zip from zenodo: https://zenodo.org/record/5730466#.YxeY4ezMLox

Due to slight numerical differences, you may not get the same flowSOM clustering exactly. You can find a vector of cluster assignments that correspond to the cells at flowsom_clusters.rds. The code has been adjusted to use these clusters.

To run the workflow:

`bash workflow.sh`

environment:
R 4.0.2

packages:
lme4 (1.1.26)
Hmisc (4.5.0)
FlowWorkspace (4.2.0)
uwot (0.1.10)
Hmisc (4.5.0)
tidyverse (1.3.0)
lmerTest (3.1.2)
FlowSOM (1.22.0)
