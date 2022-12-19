#/bin/bash

#Build the containers
for RCP in *.recipe
do
    RECIPE_NAME=$(basename $RCP)
    CONTAINER_NAME=../../tools/images/${RECIPE_NAME//'recipe'/'sif'}
    sudo singularity build ${CONTAINER_NAME} ${RECIPE_NAME} || { echo 'Container build failed, see above for details' ; exit ; } 
done