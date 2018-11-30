struct window
{
  window(){};
  window(float xmin, float xmax, float ymin, float ymax) : x1(xmin), x2(xmax), y1(ymin), y2(ymax) {};
  bool Contains(float x, float y) 
  {
    if(x>=x1 && x<=x2 && y>=y1 && y<=y2 ) 
      return true; 
    else return false; 
  };
  void Set(float xmin, float xmax, float ymin, float ymax)
  {
    x1=xmin, x2=xmax, y1=ymin, y2=ymax; 
    return;
  };
  float x1, x2, y1, y2;  
};