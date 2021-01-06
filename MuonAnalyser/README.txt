
# How to use jupyter on lxplus

## Option 1:

Just use swan.cern.ch. It mounts to eos and has all kinds of cool stuff.

## Option 2: 

Pipe jupyter through ssh. This allows us to run jupyter from /afs

Step 1:

```
ssh -L 8099:localhost:8095 <usrname>@lxplus.cern.ch
```
SSH in as normal, but send port 8095 on lxplus to your local 8099 port.

Step 2:

```
cd <dir on lxplus>
jupyter notebook --no-browser --ip=128.0.0.1 --port 8095
```

Open Jupyter and send to port.

Step 3:

on local computer, open localhost:8099 on local computer. Add in token from terminal.
