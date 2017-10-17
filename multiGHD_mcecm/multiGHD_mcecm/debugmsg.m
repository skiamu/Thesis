function debugmsg(string,variable,lines)

  %
  % Display debug message. Checks global variable DEBUG before displaying
  % a debug message.
  %
  % Code is straightfoward and also easily readable.
  %
  % -------------------------------------------------------------------
  % Author : Saket Sathe
  % Email : saket@ee.iitb.ac.in
  % Date : 9th June 2006
  % -------------------------------------------------------------------
  % 
  

  global debug
  if debug==1
    disp(string);
    disp(variable);
    if lines==1
      disp('----------------------------------------------------');
    end
  end
  
