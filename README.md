# capstone
Baruch MFE Capstone project

## Running Guide for Mac System

### Step 1: Install Command Line Tools
- Open terminal, type “xcode-select --install” in terminal (without quotes)
- A pop-up windows will appear asking you about install tools, choose install tools, wait install to finish
  
### Step 2: Install Homebrew
Run the following commands to install homebrew.
```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.bash_profile
```

Open an new terminal tab with Cmd+T (you should also close the old one), then run the following command to make sure everything works:
```
brew doctor
```

### Step 3: Install Homebrew Cask and Tap
```
brew install caskroom/cask/brew-cask
brew tap samueljohn/python
brew tap homebrew/science
```

### Step 4: Install R and its dependencies
```
brew install Caskroom/cask/xquartz
brew install R
```

### Step 5: Install Python and its modules
```
brew install python

pip install numpy
pip install scipy
pip install matplotlib
pip install pandas
```

### Step 6: Install MacTex
```
brew install Caskroom/cask/mactex
```

### Step 7: Install Github
1. Download and install the latest version of Git.
2. Tell Git your name so your commits will be properly labeled.
``` 
git config --global user.name "YOUR NAME" 
```
 
3. Tell Git the email address that will be associated with your Git commits.
```
git config --global user.email "YOUR EMAIL ADDRESS"
```
     
### Step 8: Download repository and run all
```
git clone https://github.com/rongxinyu/capstone.git
cd capstone

r capstone.r
python python-version/main.py
xelatex latex-report/hybridscheme-bss.tex
open latex-report/hybridscheme-bss.pdf
```
