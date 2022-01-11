cat $1 | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' | awk '
        BEGIN { FS=OFS="\t" }
        {for(i=2;i<=NF;i++)  
                a[$1][i]+=$i   
            a[$1][1]=i                    
        }
        END {
            for(i in a) {
                for((j=2)&&b="";j<a[i][1];j++) 
                    b=b (b==""?"":OFS)a[i][j]
                print i,b          
            }
        }' |  sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' | sed 's/ /\t/g' | awk '{for(i=4;i<=NF;i++)$i/=50}1' | sed 's/ /\t/g'