# Install npm
sudo apt-get install npm

# Install Less
sudo install less -g

# Compile prefix.less into bootstrap-iso.css will wrap .bootstrap-iso in front of everything
lessc prefix.less bootstrap-iso.css

# if you have this error /usr/bin/env: node: No such file or directory
# It's because your path needs to be updated, or if you install from a package manager you bin may be called nodejs so you just need to symlink it like so "ln -s /usr/bin/nodejs /usr/bin/node"

# TODO add script with grep to do that automatically
# You will notice that .bootstrap-iso prefixes body elements. It’s not possible for a class to prefix the body. We really want these styles to apply to just the class. We can fix this with a simple find and replace:

# Find all instance of: .bootstrap-iso body and .bootstrap-iso html
# Replace with: .bootstrap-iso

