sudo apt-get update 
sudo apt-get install -y \
    libfreetype6 fonts-dejavu fonts-dejavu-core fonts-dejavu-extra \
    fonts-liberation fonts-liberation2 \
    fontconfig \
    xvfb libxrender1 libxtst6 libxi6 libxrandr2 libxinerama1 libxfixes3 libxext6
fc-cache -f -v 
java -version 
ant -version