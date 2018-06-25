RSYNC = /usr/local/Cellar/rsync/3.1.3_1/bin/rsync
RSYNCTAGS = --archive --verbose --info=progress2
SRC_FOLDER = $(HOME)/github/dcs_summer/

REMOTE = wpq@comps0.cs.toronto.edu

# local -> comps0
synccode:
	$(RSYNC) $(RSYNCTAGS) $(SRC_FOLDER)  $(REMOTE):/u/wpq/github/dcs_summer

avoidpasswordduringssh:
	$(RSYNC) $(RSYNCTAGS) ~/.ssh/id_rsa.pub $(REMOTE):/u/wpq/.ssh
