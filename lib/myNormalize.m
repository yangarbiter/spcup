function v2 = myNormalize(v1)

v2 = v1 - min(v1) ;
v2 = v2/max(v2) ;
v2 = v2 - mean(v2) ;

end