# Troubleshooting

## Recommended steps to troubleshoot the pipeline

## Email

Check your email for an email regarding pipeline failure. You will receive an email from slurm@biowulf.nih.gov with the subject: Slurm Job_id=[#] Name=CARLISLE Failed, Run time [time], FAILED, ExitCode 1

### Review the log files

You can check logs in two ways to diagnose workflow errors:

Master SLURM Log File
Located in /path/to/results/dir/ and named slurm-[jobid].out, this file summarizes the overall workflow run. It identifies which rule failed and provides error details for any jobs executed locally via SLURM.

Individual Rule Log Files
After identifying the failed rule(s) from the master SLURM log, examine the corresponding log files in /path/to/results/dir/logs/.

### Restart the run

Each file follows this naming convention:

{rulename}.{masterjobID}.{individualruleID}.{wildcards}.{out or err}`

.out files capture standard output
.err files capture standard error messages

Once you have identified and addressed the issue, you may resume the SINCLAIR run.

Unlock the output directory, perform another dry-run, and check the status of the pipeline, then resubmit to the cluster.

```
sinclair run \
    --output /data/$USER/scRNA_test \
    -params-file assets/params.yml
```

## If a process runs out of resources

You can run sinclair with the `largemem` profile to allocate more memory and use the `largemem` slurm partition for resource-intensive processes.

```sh
sinclair run \
    -profile largemem \
    --output /data/$USER/scRNA_test \
    -params-file assets/params.yml
```

### Custom resources

You can change the resources allocated to processes by changing them in `conf/base.config`.

Here is an alternate version of the `process_high` label that allocates more resources and uses the `largemem` slurm partition:

```
    withLabel:process_high {
        cpus   = { check_max( 48                   , 'cpus'    ) }
        memory = { check_max( 500.GB * task.attempt, 'memory' ) }
        time   = { check_max( 72.h                 , 'time'    ) }
        queue = 'largemem'
        clusterOptions = ' --gres=lscratch:750 '
    }
```

## Help & Contributing

Come across a **bug**? Open an [issue](https://github.com/CCBR/SINCLAIR/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/SINCLAIR/discussions).

Want to **contribute** to this project? Check out the [contributing guidelines](../contributing.md).

**General Inquiries and Collaboration:** Please contact the CCBR Pipeliner team at [CCBR_Pipeliner@mail.nih.gov](mailto:CCBR_Pipeliner@mail.nih.gov).
