# Connecting OSC via VSCode

This tutorial will guide you through launching a SLURM job on the Ohio Supercomputer Center (OSC) and connecting to the **compute node** via VSCode Remote-SSH.


## Prerequisites
* Ensure that your SSH key is properly configured to connect to OSC.
* You will need to modify your VSCode Remote-SSH config (~/.ssh/config) accordingly.

* Update IdentityFile paths to match your SSH key locations. If you don’t have an SSH key yet, refer to this guide:  [Generating a new SSH key and adding it to the ssh-agent](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent). 
* Note that you could change `Host`, and set `IdentityFile` with the real SSH key path. you need to replace **user_name** with your **actual OSC username**.


## VSCode SSH Config Examples
You can define two hosts in your **SSH config**: one for logging into the OSC login node, and one for accessing compute nodes via a SLURM job proxy.

### For Pitzer cluster
* Add the following in ~/.ssh/config:
``` text
Host PCON0022-pitzer
    HostName pitzer.osc.edu
    User user_name 
    IdentityFile "~/.ssh/id_ed25519"

Host compute-PCON0022-pitzer
    IdentityFile "~/.ssh/id_ed25519"
    User user_name
    ProxyCommand ssh PCON0022-pitzer "~/auto_connect_pitzer.sh"
```

* Save the following `auto_connect_pitzer.sh` into your **default OSC user folder**:

``` bash
#!/bin/bash

JOB_NAME="vs_remote_pitzer" # Change to your favorite name

# Acquire job id
JOB_ID=$(squeue -h -u $USER -n $JOB_NAME -o '%A' | head -n1)

# If not exist, create a new job, modify your configuration accordingly
if [ -z "$JOB_ID" ]; then
    NODE=$(salloc --account=PCON0022 -J ${JOB_NAME} --time=05:00:00 --cpus-per-task=30 --gres=gpu:1 --no-shell 2>&1 | grep 'Nodes' | awk '{print $3}')
else
    NODE=$(squeue -h -j $JOB_ID -o '%N') # If exists, get existed job id 
fi

# Current job adress, ".ten.osc.edu" is needed.
NODE="${NODE}.ten.osc.edu"
echo $NODE
ssh -q -W ${NODE}:22 ${NODE}
```






### For Ascend cluster
* Add the following into ~/.ssh/config:
``` text
Host PCON0022-ascend
    HostName ascend.osc.edu
    User user_name 
    IdentityFile "~/.ssh/id_ed25519"

Host compute-PCON0022-ascend
    IdentityFile "~/.ssh/id_ed25519"
    User user_name
    ProxyCommand ssh PCON0022-ascend "~/auto_connect_ascend.sh"
```


* Save the following `auto_connect_ascend.sh` into your **default OSC user folder**:

``` bash
#!/bin/bash

JOB_NAME="vs_remote_ascend" # Change to your favorite name

# Acquire job id
JOB_ID=$(squeue -h -u $USER -n $JOB_NAME -o '%A' | head -n1)

# If not exist, create a new job, modify your configuration accordingly
if [ -z "$JOB_ID" ]; then
    NODE=$(salloc --account=PCON0022 -J ${JOB_NAME} --time=05:00:00 --cpus-per-task=30 --gres=gpu:1 --no-shell 2>&1 | grep 'Nodes' | awk '{print $3}')
else
    NODE=$(squeue -h -j $JOB_ID -o '%N') # If exists, get existed job id 
fi

# Current job adress, ".ten.osc.edu" is needed.
NODE="${NODE}.ten.osc.edu"
echo $NODE
ssh -q -W ${NODE}:22 ${NODE}
```


