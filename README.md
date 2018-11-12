## This is my task for my value analysis class

> ### Details  
> ### Method provided
>> #### B Plate Spline
>> #### Thin Plate Spline(TPS)
> ### The process  

The point of it is to generate a new face by converting with landmark points of other people.  
And landmark points are provided.  

### Details

Basic principle for image blending $(x^*,y^*) = f(x,y)$  
But for human it is hard to use a formula to present this process, so we try to use the 68 landmarks detected(maybe from basel model?)  

match the features point (according to the sequence)  
  
(Notice in the example provided, the background is blended and actually the landmarks' geometry is not just the same)  

### Method provided

The root method is one kind of function and after mapping the target points, i.e. the landmarks should be almost the same with provided ones, but the other pixel should make sense too.

#### B Plate Spline

> Seem that this method is not suitable for this 2D warping problem  

#### Thin Plate Spline(TPS)

Final: __I choose this principle for this problem commended by my friends__  
It's a non-rigid transformation. If we do not force the surface to pass all the provided points, i.e. with regularization, then it's something like B plate spline.  

+ [**A full Chinese reference**](https://blog.csdn.net/swimmingfish2004/article/details/7666087)  
+ [**A simple but coherent one**](https://blog.csdn.net/VictoriaW/article/details/70161180)
+ [**A tutorial exactly suits for 2D situation**](https://blog.csdn.net/kill201115/article/details/77575074)  

### The Process

Here, I will introduce the whole process of my code.  

+ Load both the images and features documents  
+ Compress the target features for Mat_L  
+ Calcu the formulations of Mat_L*Paras = Source Fea  

### Further work

It seems that I have finished most of the TPS work now at 2018.11.11. So now I have to list furthur work here.  

+ Check the __Bicubic Way__, which seems to have some waves on the face and why? - Do not know why??
+ Remove the black point on the image (the feature points creates a black hole!!!) - finished  
+ Change it to the command platte one, with users input intepolate way and the image path (or open a file searching window) and change the altered image's name with the target image index too - finished
+ Rectangle the face area and only change this area and copy this face patch back to the image - Not Needed  
+ Fasten the algorithm by use a new way for give value to image - Later

Finish the before commands before __11.12 Night__  
Later maybe:  

+ Try BPS?  
+ Finish the dlib.cpp for face feature points detection  
+ Finish the report  
