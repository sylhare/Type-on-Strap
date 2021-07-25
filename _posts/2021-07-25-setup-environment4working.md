---
published: true
tag:
  - environment
  - tools
---
## how to run my environment quickly at everywhere  

this environment only use as terminal.
# Features:
+ Multiplexer terminal with tmux, xclip for clipboard control  
+ VIM's Plugins for search/replace text, tag support cscope mapping also  
+ zsh fully support auto-complete, fuzzy finding is intergrated([fzf](https://github.com/junegunn/fzf))
+ ctags, cscope 

# 1. supported tools
+ tmux, fzf, zsh, vim (some plugin), ctags, cscope, git.

Interested VIM's Plugin:
### Vundle.vim : VIM's packet manager config is .vundle.vim
command line for install after change the config file  
``` vim +PluginInstall +q ```  

### [TagsFinder](https://github.com/AndrewRadev/tagfinder.vim) (Class, Func) you can defined how vim catch ctags result base on type, kind like this:  
```
# pri kind tag               file
  1 F C f    setOptions        SrvCfgMgr/SrvCfgTA.cpp
               class:TSrvCfgTA signature:(SPtr<TSrvParsGlobalOpt> opt)
               void TSrvCfgTA::setOptions(SPtr<TSrvParsGlobalOpt> opt)
  2 F   f    setOptions        ClntCfgMgr/ClntCfgAddr.cpp
               class:TClntCfgAddr signature:(SPtr<TClntParsGlobalOpt> opt)
               void TClntCfgAddr::setOptions(SPtr<TClntParsGlobalOpt> opt) {
  3 F   p    setOptions        ClntCfgMgr/ClntCfgAddr.h
               class:TClntCfgAddr access:public signature:(SPtr<TClntParsGlobalOpt> opt)
               void setOptions(SPtr<TClntParsGlobalOpt> opt);
```

### Files (<Space>o for open, <Space>w for save with fzf.vim)

### Find and replace quickly with followed steps:  
  refer [this](https://ttuan.github.io/2017/09/20/How-to-boost-your-vim-productivity/)  
    1.Search keyword, use /something  
    2.Type cs, replace first match, type <Esc>  
    3.Type n.n.n.n.n. To review and replace all matches
  
 
### ctags  
  ```ctags -R --c++-kinds=+p --fields=+iaS --extra=+q```
  
# 2. docker for all
  a. pull  
  	```sudo docker pull longkl/working:v1```  
  b. run  
  
  ```sudo docker  run -d  -it  --name dev_env --mount type=bind,source="/home/longkl/working"/,target="/working" -it longkl/working:v1 /bin/bash```
  