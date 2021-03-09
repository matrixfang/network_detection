#!/bin/bash
#SBATCH --job-name=liulong_zhao_sim_test              # 作业名是 
#SBATCH -N 1                 # 使用一个节点
#SBATCH --partition=compute
#SBATCH --cpus-per-task=12    # 每个进程占用12个 cpu 核心
#SBATCH -o output.out          # 将屏幕的输出结果保存到当前文件夹的 test.out
#SBATCH --mail-type=ALL # mail alert at start, end and abortion of execution
#SBATCH --mail-user=fangwenyi@pku.edu.cn #结束时发邮件给这个邮件地址

singularity exec my.sif python evolving_experiments.py
