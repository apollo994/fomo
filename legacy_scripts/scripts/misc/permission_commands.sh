chgrp -R phyloENCODE lycaenidae_benchmark # sets the group ownership
chmod -R g+rwX lycaenidae_benchmark # gives group read/write; execute only for folders and executable files
find lycaenidae_benchmark -type d -exec chmod g+s {} \; # makes new files/folders inherit the phyloENCODE group
