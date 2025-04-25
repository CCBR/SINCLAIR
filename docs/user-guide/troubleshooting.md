# Troubleshooting

Recommended steps to troubleshoot the pipeline.

## 1.1 Email

Check your email for an email regarding pipeline failure. You will receive an email from slurm@biowulf.nih.gov with the subject: Slurm Job_id=[#] Name=CARLISLE Failed, Run time [time], FAILED, ExitCode 1

## 1.2 Review the log files

Review the logs in two ways:

1. Review the master slurm file: This file will be found in the `/path/to/results/dir/` and titled `slurm-[jobid].out`. Reviewing this file will tell you what rule errored, and for any local SLURM jobs, provide error details
2. Review the individual rule log files: After reviewing the master slurm-file, review the specific rules that failed within the `/path/to/results/dir/logs/`. Each rule will include a `.err` and `.out` file, with the following formatting: `{rulename}.{masterjobID}.{individualruleID}.{wildcards from the rule}.{out or err}`

## 1.3 Restart the run

After addressing the issue, unlock the output directory, perform another dry-run and check the status of the pipeline, then resubmit to the cluster.

```
sinclair run \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --output /data/$USER/scRNA_test \
    -params-file assets/params.yml
```

## 1.4 Help & Contributing

Come across a **bug**? Open an [issue](https://github.com/CCBR/SINCLAIR/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/SINCLAIR/discussions).

Want to **contribute** to this project? Check out the [contributing guidelines](contributing.md).

**General Inquiries and Collaboration:** Please contact the CCBR Pipeliner team at [CCBR_Pipeliner@mail.nih.gov](mailto:CCBR_Pipeliner@mail.nih.gov).
