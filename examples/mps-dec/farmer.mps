NAME          FARMER   
ROWS
 N  OBJROW
 L  cons00
 G  cons01
 G  cons02
 L  cons03
 G  cons11
 G  cons12
 L  cons13
 G  cons21
 G  cons22
 L  cons23
 E  link00
 E  link01
 E  link02
 E  link10
 E  link11
 E  link12
COLUMNS
    x00       OBJROW     50             cons00     1          
    x00       cons01     3.0            link00     1  
    x01       OBJROW     76.666666      cons00     1
    x01       cons02     3.6            link01     1
    x02       OBJROW     86.666666      cons00     1
    x02       cons03    -24             link02     1
    x03       OBJROW     79.333333      cons01     1     
    x04       OBJROW     70             cons02     1           
    x05       OBJROW    -56.666666      cons01    -1        
    x06       OBJROW    -50             cons02    -1        
    x07       OBJROW    -12             cons03     1          
    x08       OBJROW    -3.333333       cons03     1          
    x10       OBJROW     50             link00    -1   
    x10       cons11     2.5            link10     1   
    x11       OBJROW     76.666666      link01    -1
    x11       cons12     3.0            link11     1
    x12       OBJROW     86.666666      link02    -1
    x12       cons13    -20             link12     1
    x13       OBJROW     79.333333      cons11     1     
    x14       OBJROW     70             cons12     1           
    x15       OBJROW    -56.666666      cons11    -1        
    x16       OBJROW    -50             cons12    -1        
    x17       OBJROW    -12             cons13     1          
    x18       OBJROW    -3.333333       cons13     1          
    x20       OBJROW     50             link10    -1       
    x20       cons21     2.0            
    x21       OBJROW     76.666666      link11    -1
    x21       cons22     2.4            
    x22       OBJROW     86.666666      link12    -1
    x22       cons23    -16             
    x23       OBJROW     79.333333      cons21     1     
    x24       OBJROW     70             cons22     1           
    x25       OBJROW    -56.666666      cons21    -1        
    x26       OBJROW    -50             cons22    -1        
    x27       OBJROW    -12             cons23     1          
    x28       OBJROW    -3.333333       cons23     1          
RHS
    RHS1      cons00     500            cons01     200       
    RHS1      cons02     240            cons11     200       
    RHS1      cons12     240            cons21     200       
    RHS1      cons22     240  
BOUNDS
 UP BOUND     x00        1e+30
 UP BOUND     x01        1e+30
 UP BOUND     x02        1e+30
 UP BOUND     x07        6000       
 UP BOUND     x10        1e+30
 UP BOUND     x11        1e+30
 UP BOUND     x12        1e+30
 UP BOUND     x17        6000       
 UP BOUND     x20        1e+30
 UP BOUND     x21        1e+30
 UP BOUND     x22        1e+30
 UP BOUND     x27        6000       
ENDATA
