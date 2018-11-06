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
