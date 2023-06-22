# dotfiles
git clone https://github.com/jeksterslab/dotfiles.git
cd dotfiles
(cd .setup/bash && make) || echo Building .bashrc failed.
# (cd .setup/git && make) || echo Building .gitconfig failed.
# (cd .setup/ssh && make) || echo Building ssh keys failed.
(cd .setup/tmux && make) || echo Building .tmux.conf failed.
(cd .setup/vim && make) || echo Building .vimrc failed.
# (cd .setup/xterm && make) || echo Building .Xresources failed.
cd ..
rm -rf dotfiles
