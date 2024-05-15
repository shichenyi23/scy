<!--
 * @Author: shichenyi shichenyi23@163.com
 * @Date: 2024-05-15 16:16:58
 * @LastEditors: shichenyi shichenyi23@163.com
 * @LastEditTime: 2024-05-15 17:14:00
 * @FilePath: \hw1\readme\README.md
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
-->
# algebra
硬件技术团队编程基础作业
## 预修要求｜Requirements
修读过《C程序设计基础》、《线性代数》以及X·Lab硬件技术团队第一次内训或其对应的高阶课程。
## 说明｜Explainations
本题目的难度对于初学者而言相对较高，主要考察了基础的数学能力、通过编程解决问题的能力以及工程管理、CMake、git等综合能力。该作业的最终得分仅作参考，同学们可根据自己的能力来决定实现哪些函数。
## 题目背景｜Background
《线性代数》作为浙江大学工科多数专业必修的数学基础课程，对于其掌握是至关重要的，后续各大专业的专业课程也都离不开线性代数。然而，在后续的专业课程学习中，我们往往只需要计算一些矩阵的数值解，这个过程如果用手去计算的话是十分痛苦的。秉承着“我都学编程了就不要自己做一些无意义的事情”的原则，我们决定实现一个线性代数计算库，来辅助我们进行运算。
> 当然，如MATLAB、Python等高级编程语言已经可以做到这类事情，且做得更好，但这并不妨碍我们通过这样一种方式来锻炼自己的C语言编程能力。
## 题目介绍｜Introduction
本仓库给出了我们在内训中提到的工程模板，同学们要做的任务如下：
1. 自学git，注册GitHub账号，将本仓库在自己的GitHub账户下Fork一份（注意是Fork，禁止直接clone本仓库到本地，否则你将无法完成后续提交），并按照`yourname_hw1`的格式更改仓库名称（在仓库中的Settings处可修改，记得不要用中文，仓库权限为public，如涉及到隐私保护，可设为private，但要将`tanninrachel@yinlin.wiki`这个账户设置为协作者）。
2. 将你的仓库clone到本地。
3. 按照内训所讲的工程模板补充所缺的文件夹。
4. 根据`inc/algebra.h`中的注释和预定义，在`src/algebra.c`中实现对应的函数。
5. 根据内训所讲，自行编写`CMakeLists.txt`文件，使你的工程能够在本地成功编译运行。
6. 自学Markdown，修改`README.md`文件，需要包含你的实现思路（大致描述即可）以及本地运行截图。
7. 将你的修改提交到远程仓库，并将仓库链接提交（提交方式待定）。
## 思路参考｜Thinking
见`doc`文件夹。
## 交互格式｜Format
在本题目中，`main.c`文件已给出，不需要同学们自己实现，也请不要更改这个文件，否则可能出现判题错误。
### 输入格式
本题目采用帧判定的思路进行，每一帧的第一行指令代码，`+`、`-`、`*`、`.`、`t`、`d`、`i`、`r`、`j`分别测试`add_matrix`、`sub_matrix`、`mul_matrix`、`scale_matrix`、`transpose_matrix`、`det_matrix`、`inv_matrix`、`rank_matrix`、`trace_matrix`函数，`q`表示退出。

接下来的一行输入矩阵 $\mathbf{A}$ 的行数 $m$ 和列数 $n$ ，在接下来的 $m$ 行中输入 $n$ 个双精度浮点数，以空格分开。
对于二元运算函数的测试，需要再按照上述过程输入矩阵 $\textbf{B}$ 。

可能的一次运行输入如下：
```
+
2 2
1.1 1.3
2.4 3.7
2 2
3.1 4.3
5.1 7.1
+
2 2
1.1 1.2
2.4 3.5
2 3
1 2 2.1
3 2 3.3
q
```
### 输出格式
在每一帧中，依次根据输入的指令代码运行对应的函数，并给出函数的输出与标准值比对。上述输入的正确输出如下：
```
4.20    5.60    
7.50    10.80
Error: Matrix a and b must have the same rows and cols.
```
## 评分标准｜Standard
* 成功运行：+25分
* `add_matrix`、`sub_matrix`、`mul_matrix`、`scale_matrix`、`transpose_matrix`、`trace_matrix`功能正常每个+5分
* `det_matrix`、`inv_matrix`、`rank_matrix`功能正常每个+15分
> 为了保证题目难度，每个函数具体的评分标准不予公布。
## 提交方式｜Submit
将你完成任务的GitHub仓库链接通过飞书发送给X-Check机器人，机器人会自动进行打分并将分数返回。参考如下视频：

[https://tannin-1316822731.cos.ap-nanjing.myqcloud.com/2024-03-29-827_1711724348.mp4](https://tannin-1316822731.cos.ap-nanjing.myqcloud.com/2024-03-29-827_1711724348.mp4)
## 实现思路及运行截图
# `add_matrix`
![alt text](<屏幕截图 2024-05-15 163526.png>)
# `sub_matrix`
![alt text](<屏幕截图 2024-05-15 165715.png>)
# `mul_matrix`
![alt text](<屏幕截图 2024-05-15 165715-1.png>)
# `scale_matrix`
![alt text](<屏幕截图 2024-05-15 165715-2.png>)
# `transpose_matrix`
![alt text](<屏幕截图 2024-05-15 165741.png>)
# `trace_matrix`
![alt text](<屏幕截图 2024-05-15 165826.png>)
# `det_matrix`
行列式运用的是按第一行展开，至于余子式的计算，先将其中一行一列去掉后，重新调用这个函数进行计算，直到只剩一行一列，行列式的值为其本身，要用到递归思想。
![alt text](<屏幕截图 2024-05-15 170215.png>)
# `inv_matrix`
算矩阵的逆，先调用det_matrix函数，判断有没有逆，若不等于0，先算伴随矩阵，再除以行列式
![alt text](<屏幕截图 2024-05-15 170222.png>)
# `rank_matrix`
运用高斯消元法，注意交换0所在行
![alt text](<屏幕截图 2024-05-15 171257.png>)
## 截止日期｜Deadline
待定