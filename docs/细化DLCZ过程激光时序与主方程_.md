# **基于光纤腔的2到3个中性原子的DLCZ过程的详细分析**

**1\. 导言**

Duan-Lukin-Cirac-Zoller (DLCZ) 协议在实现远距离量子通信和分布式量子计算方面发挥着至关重要的作用。该协议利用原子系综作为量子存储器，通过光子作为信息载体，实现了量子纠缠的远距离传输。中性原子因其相干时间长而成为理想的量子存储介质，而光纤腔则通过增强光与物质的相互作用，显著提高了光子发射和吸收的效率，这对于高效地实现 DLCZ 协议至关重要 1。光纤腔的小模式体积和高精细度能够极大地增强单个原子与腔内光场之间的耦合强度 2，从而实现更高效的量子操作，这对于研究和应用 DLCZ 协议至关重要。本报告旨在对基于光纤腔的 2 到 3 个中性原子的 DLCZ 过程的六个步骤进行详细分析，重点关注每个步骤中涉及的激光控制、潜在的噪声来源以及描述系统演化的主方程。通过深入探讨近期文献中的实验进展和理论描述，本报告将为从事或计划从事基于光纤腔的少原子 DLCZ 研究的实验量子光学专家提供全面的信息。

**2\. DLCZ过程的六个步骤回顾**

为了确保后续分析的连贯性，本节将简要回顾 DLCZ 协议通常包含的六个基本步骤。这些步骤为后续详细讨论奠定了基础。虽然用户未提供之前的报告，但根据标准的 DLCZ 协议，我们可以假设这六个步骤包括：(1) 将本地和远程的少数中性原子（2-3 个，一侧至少有一个）捕获并冷却到光纤腔的模式体积内；(2) 将原子初始化到特定的量子态；(3) 使用弱写入激光脉冲激发原子，并在光纤腔中发射一个信号光子；(4) 将原子激发（通常是自旋波或亚稳态）存储一段时间；(5) 使用读出激光脉冲读取原子激发，并在光纤腔中发射一个闲置光子；(6) 通过对来自不同原子（或原子系综）的信号光子和闲置光子进行干涉和检测，实现远距离原子之间的纠缠。后续章节将根据最新的研究进展，对这些步骤进行更详尽的阐述和完善。

**3\. 基于近期文献的DLCZ过程各步骤的详细分析**

本节将深入分析 DLCZ 过程的六个步骤，并结合提供的研究摘要，详细阐述每个步骤的实验流程、理论描述以及噪声考虑。

**3.1 步骤1：在中空光纤腔中捕获和冷却中性原子**

* **过程描述与激光时序：** 该初始步骤涉及将少量中性原子（2-4 个）捕获并冷却到光纤腔的模式体积内。高效的捕获和冷却对于实现强的原子-腔相互作用和长的相干时间至关重要。通常，原子首先通过磁光阱 (MOT) 从低压蒸汽中捕获 2，MOT 利用多个反向传播的激光束（频率略低于原子共振频率）和一个四极磁场。随后，原子通常被转移到由强聚焦激光束（光镊）或驻波形成的光学偶极阱 (ODT) 中 1。Snippet 1 讨论了使用沿腔轴或垂直于腔轴的光学偶极阱将中性原子限制在腔镜之间。Snippet 2 描述了使用光阱将冷原子转移到光纤腔中。Snippet 5 强调了光镊在捕获单个原子方面的应用以及激光强度噪声在引起参量加热和限制原子寿命方面的关键作用。ODT 激光的波长通常远失谐于任何原子跃迁，以最大限度地减少散射。强度决定了阱深，加载持续时间取决于转移效率。为了进一步降低 ODT 中原子的温度，通常会采用额外的冷却阶段，例如光学粘滞冷却或亚多普勒冷却技术 2。Snippet 2 提到在 MOT 捕获后应用光学粘滞冷却将原子进一步冷却至约 5 μK。Snippet 6 描述了使用沿重力方向的推束来减少单个原子加载到腔模式的时间并缩小其速度分布。  
* **哈密顿量 (LaTeX):** 此步骤的哈密顿量主要描述原子与捕获激光场的相互作用，可以用势能项来建模。对于 ODT，势能与激光束的强度成正比。  
  * H=∑i=1N​(2mpi2​​+Vtrap​(ri​))，其中 N 是原子数，$ \\mathbf{p}*i $ 是第 i 个原子的动量，$ m $ 是原子质量，$ V*{trap}(\\mathbf{r}*i) $ 是第 i 个原子位置处的捕获势能。$ V*{trap} $ 的形式取决于阱的几何形状（例如，对于光镊是高斯型）。  
* **Lindblad 算符 (Lk​) 及对应的速率 (γk​) (LaTeX):**  
  * **激光强度噪声引起的加热：** 正如 Snippets 5 所讨论的，在两倍于阱频率的频率下，捕获激光的强度波动会导致原子的参量加热，增加其动能并可能导致其逃离阱。这可以用增加原子动量或振动能量的 Lindblad 算符来建模。速率与相关频率下激光噪声的功率谱密度成正比。  
  * **与背景气体的碰撞：** Snippets 1 表明，与真空室中残留气体分子的碰撞会导致原子损失。此过程可以用以速率与背景气体密度和碰撞截面成正比的方式从系统中移除原子的 Lindblad 算符来描述。Snippet 7 专门关注背景碰撞引起的损失速率常数。  
  * **阱激光的自发辐射：** 如果捕获激光的失谐不够大，则可能发生自发辐射，导致原子加热或损失。这可以用标准的原子衰变 Lindblad 算符建模，速率取决于激光的失谐和强度。  
* **主方程 (LaTeX):**  
  * dtdρ​=−ℏi​\[H,ρ\]+∑i=1N​γheat,i​(Aheat,i​ρAheat,i†​−21​{Aheat,i†​Aheat,i​,ρ})+γloss​(Lloss​ρLloss†​−21​{Lloss†​Lloss​,ρ})+∑i=1N​γspont,trap,i​(σtrap,i−​ρσtrap,i+​−21​{σtrap,i+​σtrap,i−​,ρ})  
  * 其中，$ A\_{heat, i} $ 是描述第 i 个原子加热的算符，$ \\gamma\_{heat, i} $ 是加热速率，$ L\_{loss} $ 是原子损失算符，$ \\gamma\_{loss} $ 是由于背景气体碰撞造成的损失速率，$ \\sigma^-*{trap, i} $ 和 $ \\sigma^+*{trap, i} $ 是捕获激光驱动的原子跃迁的降低和升高算符，$ \\gamma\_{spont, trap, i} $ 是由于捕获激光造成的自发辐射速率。  
* **主要噪声：** 激光强度噪声（参量加热）、与背景气体的碰撞（原子损失）、阱激光的自发辐射（加热/损失）。

稳定的捕获和低温是 DLCZ 协议后续步骤的先决条件。捕获激光中的噪声是限制捕获原子寿命的一个重要因素。由激光强度波动驱动的参量加热会增加捕获原子的能量，使其更容易逃脱捕获势阱。这直接减少了可用于 DLCZ 过程的原子数量，并限制了实验的持续时间。Snippets 5 中关于降低激光强度噪声及其对光镊中原子寿命影响的研究强调了这一点的重要性。

**3.2 步骤2：原子态的初始化**

* **过程描述与激光时序：** 一旦将所需数量的原子（2-3 个）捕获并冷却到光纤腔的模式体积内，就需要将它们初始化到明确定义的量子态，通常是与 DLCZ 协议相关的电子基态。这通常涉及应用一个或多个调谐到特定原子跃迁的激光束，以将布居转移到所需的基态。激光通常是圆偏振的，并在弱磁场存在下应用以打破简并。泵浦激光的持续时间和强度经过优化，以实现高保真度的初始化。  
* **哈密顿量 (LaTeX):** 光学泵浦的哈密顿量描述了原子与泵浦激光场的相互作用，诱导不同原子能级之间的跃迁。  
  * H=∑i=1N​ℏΩpump,i​(σpump,i+​e−iωpump​t+σpump,i−​eiωpump​t)，其中 $ \\Omega\_{pump, i} $ 是第 i 个原子的泵浦激光的拉比频率，$ \\sigma^+\_{pump, i} $ 和 $ \\sigma^-*{pump, i} $ 是相关跃迁的升高和降低算符，$ \\omega*{pump} $ 是激光频率。  
* **Lindblad 算符 (Lk​) 及对应的速率 (γk​) (LaTeX):**  
  * **泵浦过程中的自发辐射：** 当原子被泵浦激光激发时，它们会自发辐射回到较低的能级。虽然理想情况下它们会衰变到目标基态，但也有可能衰变到其他态，从而降低泵浦过程的效率。这可以用标准的原子自发辐射 Lindblad 算符建模，速率对应于跃迁概率。  
  * **离共振散射：** 泵浦激光可能与其他原子能级存在微弱的离共振相互作用，导致不希望的激发或加热。这些可以用额外的 Lindblad 算符建模。  
* **主方程 (LaTeX):**  
  * dtdρ​=−ℏi​\[H,ρ\]+∑i=1N​∑j​γspont,pump,ij​(σpump,ij−​ρσpump,ij+​−21​{σpump,ij+​σpump,ij−​,ρ})  
  * 其中，$ \\sigma^-*{pump, ij} $ 和 $ \\sigma^+*{pump, ij} $ 是第 i 个原子的光学泵浦中涉及的第 j 个跃迁的降低和升高算符，$ \\gamma\_{spont, pump, ij} $ 是相应的自发辐射速率。  
* **主要噪声：** 不完善的光学泵浦（布居停留在不需要的状态）、自发辐射到非目标态、离共振散射。

高保真度的初始化对于后续的纠缠生成步骤至关重要。激发态或其他基态中的任何残留布居都可能导致 DLCZ 过程中的错误。DLCZ 协议依赖于从明确定义的初始态开始的特定拉曼跃迁。如果原子没有被正确初始化，写入脉冲可能会诱导不希望的跃迁，导致噪声光子的产生或在不正确的态中创建原子激发，从而降低所生成纠缠的保真度。

**3.3 步骤3：写入原子激发并向光纤腔发射信号光子**

* **过程描述与激光时序：** 将一个弱写入激光脉冲施加到捕获并初始化的原子上。该脉冲被调谐以诱导拉曼跃迁，从而在一个原子中产生激发（或多个原子中激发的叠加），并同时将一个信号光子发射到光纤腔的模式中。脉冲的弱性确保了产生超过一个激发-光子对的概率很低。写入激光脉冲的频率失谐于中间激发态，使其与腔模式共同驱动从初始基态到目标态（例如，另一个基态或亚稳态）的拉曼跃迁。脉冲持续时间应短于中间激发态的寿命，以最大限度地减少自发辐射。强度保持较低，以确保单激发状态。Snippet 8 提到了一个产生虚拟能级的强写入/读取脉冲，这在 DLCZ 中是一种常见技术。Snippets 9 描述了一个诱导拉曼跃迁并发射光子的弱写入激光。Snippet 11 提到了用于时间多路复用的写入脉冲序列。  
* **哈密顿量 (LaTeX):** 写入过程的哈密顿量包括原子与写入激光场和光纤腔的量子化场的相互作用。这可以用拉曼过程产生的有效相互作用哈密顿量来描述。  
  * Hint​=∑i=1N​ℏ(gi​a†σeg,i−​+Ωwrite,i​σse,i+​)+h.c. （这是一个简化的形式，侧重于相关的跃迁；更详细的哈密顿量将包括失谐和完整的能级结构）。其中，$ a^{\\dagger} $ 是腔模式中光子的产生算符，$ \\sigma^-*{eg, i} $ 是原子 i 从激发态 $ |e\\rangle $ 到基态 $ |g\\rangle $ 的降低算符（与腔耦合），$ \\Omega*{write, i} $ 是原子 i 从初始态 $ |s\\rangle $ 到激发态 $ |e\\rangle $ 的跃迁的写入激光的拉比频率，$ \\sigma^+\_{se, i} $ 是相应的升高算符。$ g\_i $ 是原子 i 与腔的耦合强度。  
* **Lindblad 算符 (Lk​) 及对应的速率 (γk​) (LaTeX):**  
  * **中间激发态的自发辐射：** 被写入激光激发到中间能级的原子会自发衰变到其他态，导致相干性损失和潜在的错误。这可以用标准的原子自发辐射 Lindblad 算符建模，速率对应于衰变通道。  
  * **腔衰变：** 发射到腔中的光子会由于光纤腔镜的不完美或光纤内的吸收而损失。这由腔衰变 Lindblad 算符描述，速率为 $ \\kappa $。  
  * **原子退相干：** 在写入过程中，原子态会由于与环境的相互作用而发生退相干。这可以用适当的 Lindblad 算符（例如，退相位）建模。  
* **主方程 (LaTeX):**  
  * dtdρ​=−ℏi​\[Hint​,ρ\]+∑i=1N​∑j​γspont,write,ij​(σwrite,ij−​ρσwrite,ij+​−21​{σwrite,ij+​σwrite,ij−​,ρ})+κ(aρa†−21​{a†a,ρ})  
  * 其中，$ \\gamma\_{spont, write, ij} $ 是第 i 个原子的第 j 个衰变通道中中间激发态的自发辐射速率。  
* **主要噪声：** 中间激发态的自发辐射、腔衰变、脉冲期间的原子退相干。

以高保真度生成单个激发和单个光子的效率至关重要。写入过程中发生的自发辐射是需要最小化的主要噪声源。自发辐射会导致原子在拉曼跃迁完成之前衰变到不需要的状态，或者导致发射的光子具有与所需信号光子不同的频率或模式。这降低了成功预示的概率，并引入了后续纠缠生成中的错误。正如 Snippet 8 中提到的，使用远失谐拉曼跃迁有助于抑制荧光噪声。

**3.4 步骤4：原子激发（自旋波）的存储**

* **过程描述与激光时序：** 在成功检测到信号光子（预示原子激发的产生）后，这种激发（通常是自旋波或亚稳态）需要存储一段时间。存储时间受各种退相干机制的限制。理想情况下，在存储期间不施加激光。然而，在某些实现中，可能会使用弱保持激光来防止激发的扩散或在存储期间操纵原子态。这些激光的参数将取决于具体的实现。  
* **哈密顿量 (LaTeX):** 存储期间的哈密顿量理想情况下描述了原子激发的自由演化。然而，在现实中，由于残留磁场或与捕获激光的相互作用，可能存在小的能级移动。  
  * H=∑i=1N​ℏδi​σstore,iz​，其中 $ \\delta\_i $ 是原子 i 的存储态的一个小失谐，$ \\sigma^z\_{store, i} $ 是存储态的泡利 Z 算符。  
* **Lindblad 算符 (Lk​) 及对应的速率 (γk​) (LaTeX):**  
  * **存储态的退相干：** 存储的原子激发会由于各种因素而发生退相干，包括与环境的相互作用、磁场波动以及自发跃迁到其他长寿命态。这可以用导致相干性衰减（例如，相位翻转或幅度阻尼）的 Lindblad 算符建模，速率为 $ \\gamma\_{dephasing} $。Snippets 12 讨论了捕获原子中的退相干机制，例如由捕获噪声引起的声子跳跃诱导退相干 (PJID) 和由微分光移 (DLS) 引起的退相干，这些可能与此相关。Snippets 13 提到了光子逃逸和原子衰变引起的退相干。Snippet 14 指出了弹性碰撞和非弹性碰撞引起的退相干。  
  * **原子激发的损失：** 在某些情况下，存储的激发可能会由于跃迁到其他态而丢失。这可以用幅度阻尼 Lindblad 算符建模，速率为 $ \\gamma\_{loss, store} $。  
* **主方程 (LaTeX):**  
  * dtdρ​=−ℏi​\[H,ρ\]+∑i=1N​γdephasing,i​(σstore,iz​ρσstore,iz​−ρ)+∑i=1N​γloss,store,i​(σstore,i−​σstore,i+​ρ−ρσstore,i−​σstore,i+​) （这是一个简化的例子；确切的形式取决于存储激发的具体类型和主要的退相干机制）。  
* **主要噪声：** 存储自旋波的退相干（退相位、幅度阻尼）、原子激发的损失。

存储原子激发的相干时间是 DLCZ 量子存储器性能的关键参数。更长的相干时间允许在更远的距离上进行纠缠，并进行更复杂的量子操作。退相干会降低存储在原子激发中的量子信息的质量，从而降低后续操作（如纠缠交换或隐形传态）的保真度。存储时间限制了使用此协议可靠传输量子信息的距离。环境噪声（如磁场波动或捕获势的不稳定性）会对存储原子态的退相干产生重大影响。这些波动会导致原子能级的能量移动，从而导致代表自旋波的叠加态的相位相干性损失。最大限度地减少这些环境噪声对于延长相干时间至关重要。

**3.5 步骤5：读取原子激发并向光纤腔发射闲置光子**

* **过程描述与激光时序：** 为了检索存储的量子信息，需要向原子施加一个读出激光脉冲。该脉冲激发一个拉曼跃迁，将原子激发重新转换为一个光子，该光子作为闲置光子发射到光纤腔模式中。读出激光脉冲的频率被选择为与一个跃迁共振，该跃迁与腔模式共同作用，可以将存储的原子激发转换为光子。强度和持续时间经过优化，以实现高效和相干的检索。Snippet 8 提到了一个强写入/读取脉冲。Snippet 9 描述了一个读取激光，用于将自旋波激发检索为闲置光子。  
* **哈密顿量 (LaTeX):** 读取过程的哈密顿量描述了原子与读取激光和腔模式的相互作用，同样通过有效的拉曼相互作用。  
  * Hint​=∑i=1N​ℏ(gi​a†σgs,i−​+Ωread,i​σes,i+​)+h.c. （同样，这是一个简化的形式）。其中，$ \\sigma^-*{gs, i} $ 是原子 i 从激发态 $ |e\\rangle $ 到存储态 $ |s\\rangle $ 的降低算符（在读取期间与腔耦合），$ \\Omega*{read, i} $ 是原子 i 从基态 $ |g\\rangle $ 到激发态 $ |e\\rangle $ 的跃迁的读取激光的拉比频率，$ \\sigma^+\_{es, i} $ 是相应的升高算符。  
* **Lindblad 算符 (Lk​) 及对应的速率 (γk​) (LaTeX):**  
  * **读取过程中的自发辐射：** 与写入过程类似，在读取脉冲期间，中间激发态的自发辐射会引入噪声并降低检索效率。  
  * **腔衰变：** 发射的闲置光子可能会从腔中丢失。  
  * **原子退相干：** 在读取脉冲期间的任何剩余退相干也会影响检索光子的保真度。  
* **主方程 (LaTeX):**  
  * dtdρ​=−ℏi​\[Hint​,ρ\]+∑i=1N​∑j​γspont,read,ij​(σread,ij−​ρσread,ij+​−21​{σread,ij+​σread,ij−​,ρ})+κ(aρa†−21​{a†a,ρ})  
  * 其中，$ \\gamma\_{spont, read, ij} $ 是读取过程中中间激发态的自发辐射速率。  
* **主要噪声：** 读取过程中的自发辐射、腔衰变、脉冲期间的原子退相干。

高效且相干地将存储的激发检索为闲置光子对于 DLCZ 协议的整体成功至关重要。此步骤中的损失和噪声会显着降低所生成纠缠的保真度。读取脉冲中的任何缺陷或读取过程中发生的自发辐射事件都可能导致原子激发被转换为具有不正确属性（例如，错误的偏振、频率或时间模式）的光子。这直接影响了在远距离原子存储器之间建立的纠缠质量。如果在读取过程完成之前原子发生自发衰变，则存储的量子信息就会丢失。这强调了使用适当的原子能级方案和激光参数以最大限度地减少自发辐射的重要性。

**3.6 步骤6：通过光子干涉和探测实现远距离原子间的纠缠**

* **过程描述与激光时序：** 从一个原子系综（或几个原子）发射的信号光子和从另一个远距离原子系综（或几个原子）发射的闲置光子被发送到中央测量站。然后，这些光子在分束器上发生干涉，并使用探测器记录光子的到达。探测事件预示着两个空间分离的原子激发之间的纠缠。可能需要重复此过程以达到所需的纠缠概率。此步骤主要涉及发射光子的操作和检测，没有直接的激光脉冲施加到原子上。  
* **哈密顿量 (LaTeX) 和 Lindblad 算符 (Lk​) 及对应的速率 (γk​) (LaTeX):** 这些与原子在此步骤中的状态没有直接关系，而是描述了光子在光路中的演化和探测过程。  
* **主方程 (LaTeX):** 不直接适用于原子状态。  
* **主要噪声：**  
  * **传输通道（例如，光纤）中的光子损失：** 当光子通过光纤传输时，会被吸收或散射，导致损失 15。Snippet 16 讨论了在 12 公里光纤上生成纠缠，强调了光纤损耗带来的挑战。  
  * **分束器处的不完美干涉：** 为了成功纠缠，两个光子必须是不可区分的，并在分束器处完美干涉。它们的时间模式、频率或偏振的任何差异都会降低干涉的可见度。  
  * **探测器效率低下和暗计数：** 现实世界中的光子探测器的效率低于 100%，这意味着一些光子可能无法被检测到。探测器还具有非零的暗计数率，即使没有光子存在，也会记录到探测事件，从而导致误报。Snippet 17 提到了单光子计数模块的暗计数率。  
  * **光路中的相位波动：** 光纤或其他光学元件折射率的波动会导致两个光子之间的相位差，从而影响它们的干涉。  
* **来自 Snippets 的见解：** Snippets 8 都与使用原子和光子生成和分配纠缠有关，通常涉及腔。它们强调了该过程的概率性以及影响纠缠生成保真度和速率的各种因素。Snippet 8 提到了低无条件噪声和违反柯西-施瓦茨不等式是高保真度生成的指标。Snippet 16 强调了纠缠生成速率需要快于存储器退相干速率。Snippets 21 讨论了纳米光纤腔在快速纠缠生成方面的潜力。Snippets 22 展示了使用光纤腔的高保真度量子比特读出和纠缠生成。

使用 DLCZ 协议纠缠远距离原子本质上是概率性的，因为它需要成功检测传输和干涉后的光子。纠缠的速率和保真度受到光通道和探测器损失和噪声的显着影响。纠缠的预示依赖于特定的光子检测事件。连接远距离原子到中央测量站的光纤中的损耗降低了这些光子到达探测器的概率。不完美的干涉或有噪声的探测器会导致对原子状态的错误推断，从而限制了可以建立高保真度纠缠的速率。

**4\. 精细化DLCZ过程的总结表**

本节将提供一个表格，总结每个步骤的详细分析，包括精细化的过程描述、激光脉冲序列和参数、哈密顿量和 Lindblad 算符（以 LaTeX 格式）以及主要噪声来源。

| 步骤 | 过程描述 | 激光脉冲序列和参数 | 哈密顿量 (LaTeX) | Lindblad 算符及速率 (LaTeX) | 主要噪声 |
| :---- | :---- | :---- | :---- | :---- | :---- |
| 1 | 捕获和冷却 | MOT 激光、ODT 激光（波长、强度、持续时间）、冷却激光（参数）、推束（可选） | H=∑i=1N​(2mpi2​​+Vtrap​(ri​)) | γheat,i​(Aheat,i​ρAheat,i†​−21​{Aheat,i†​Aheat,i​,ρ}), γloss​(Lloss​ρLloss†​−21​{Lloss†​Lloss​,ρ}), γspont,trap,i​(σtrap,i−​ρσtrap,i+​−21​{σtrap,i+​σtrap,i−​,ρ}) | 激光强度噪声、背景气体碰撞、阱激光的自发辐射 |
| 2 | 初始化 | 光学泵浦激光（波长、偏振、强度、持续时间）、弱磁场 | H=∑i=1N​ℏΩpump,i​(σpump,i+​e−iωpump​t+σpump,i−​eiωpump​t) | γspont,pump,ij​(σpump,ij−​ρσpump,ij+​−21​{σpump,ij+​σpump,ij−​,ρ}) | 不完善的光学泵浦、自发辐射到非目标态、离共振散射 |
| 3 | 写入和信号光子发射 | 弱写入激光脉冲（频率、强度、持续时间、失谐） | Hint​=∑i=1N​ℏ(gi​a†σeg,i−​+Ωwrite,i​σse,i+​)+h.c. | γspont,write,ij​(σwrite,ij−​ρσwrite,ij+​−21​{σwrite,ij+​σwrite,ij−​,ρ}), κ(aρa†−21​{a†a,ρ}) | 中间态的自发辐射、腔衰变、原子退相干 |
| 4 | 激发存储 | （理想情况下没有激光）或保持激光（参数） | H=∑i=1N​ℏδi​σstore,iz​ | γdephasing,i​(σstore,iz​ρσstore,iz​−ρ), γloss,store,i​(σstore,i−​σstore,i+​ρ−ρσstore,i−​σstore,i+​) | 存储态的退相干、原子激发的损失 |
| 5 | 读取和闲置光子发射 | 读取激光脉冲（频率、强度、持续时间） | Hint​=∑i=1N​ℏ(gi​a†σgs,i−​+Ωread,i​σes,i+​)+h.c. | γspont,read,ij​(σread,ij−​ρσread,ij+​−21​{σread,ij+​σread,ij−​,ρ}), κ(aρa†−21​{a†a,ρ}) | 读取过程中的自发辐射、腔衰变、原子退相干 |
| 6 | 通过光子干涉实现纠缠 | 信号光子和闲置光子干涉后的探测 | 不直接适用于原子态 | 不直接适用于原子态 | 通道中的光子损失、不完美的干涉、探测器效率低下和暗计数、相位波动 |

**5\. 结论**

本报告详细分析了基于光纤腔的 2 到 3 个中性原子的 DLCZ 过程的六个步骤，重点介绍了关键的激光参数、描述底层量子动力学的主方程以及每个步骤的主要噪声来源。近期的研究进展表明，在使用该协议实现高保真度纠缠生成和分发方面取得了显著进展。未来的研究方向可能包括探索新的降噪技术 30，或提高基于 DLCZ 的量子网络的效率和可扩展性 21。关于腔暗模式和集体效应的研究 32 也可能为提高性能提供途径。

#### **Works cited**

1. Cavity quantum electrodynamics \- The domain cambridgecore.org is registered by NetNames, accessed March 22, 2025, [https://core-cms.cambridgecore.org/core/services/aop-cambridge-core/content/view/923744B83152A17880D819D19B34102E/9780511762314c17\_p280-295\_CBO.pdf/cavity\_quantum\_electrodynamics.pdf](https://core-cms.cambridgecore.org/core/services/aop-cambridge-core/content/view/923744B83152A17880D819D19B34102E/9780511762314c17_p280-295_CBO.pdf/cavity_quantum_electrodynamics.pdf)  
2. Cavity Quantum Electrodynamics (QED), accessed March 22, 2025, [https://www.mpq.mpg.de/4939125/qed](https://www.mpq.mpg.de/4939125/qed)  
3. Deterministic Generation of Single Photons from One Atom Trapped in a Cavity | Request PDF \- ResearchGate, accessed March 22, 2025, [https://www.researchgate.net/publication/8682388\_Deterministic\_Generation\_of\_Single\_Photons\_from\_One\_Atom\_Trapped\_in\_a\_Cavity](https://www.researchgate.net/publication/8682388_Deterministic_Generation_of_Single_Photons_from_One_Atom_Trapped_in_a_Cavity)  
4. Introduction of Cavity QED \- Nanofiber Quantum Technologies, accessed March 22, 2025, [https://www.nano-qt.com/resources/20241101001](https://www.nano-qt.com/resources/20241101001)  
5. (PDF) Reduction of laser intensity noise over 1 MHz band for single ..., accessed March 22, 2025, [https://www.researchgate.net/publication/345471909\_Reduction\_of\_laser\_intensity\_noise\_over\_1\_MHz\_band\_for\_single\_atom\_trapping](https://www.researchgate.net/publication/345471909_Reduction_of_laser_intensity_noise_over_1_MHz_band_for_single_atom_trapping)  
6. \[2403.03019\] Pushing single atoms near an optical cavity \- arXiv, accessed March 22, 2025, [https://arxiv.org/abs/2403.03019](https://arxiv.org/abs/2403.03019)  
7. Study of background gas collisions in atomic traps \- Madison Group, accessed March 22, 2025, [https://qdglab.physics.ubc.ca/files/2022/05/thesis\_phd\_VanDongen.pdf](https://qdglab.physics.ubc.ca/files/2022/05/thesis_phd_VanDongen.pdf)  
8. A broadband DLCZ quantum memory in room-temperature atoms, accessed March 22, 2025, [https://d-nb.info/1170134726/34](https://d-nb.info/1170134726/34)  
9. Experimental entanglement of 25 individually accessible atomic ..., accessed March 22, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC5930417/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5930417/)  
10. Experimental realization of a multiplexed quantum memory with 225 individually accessible memory cells \- PMC, accessed March 22, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC5424256/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5424256/)  
11. Cavity-enhanced and temporally multiplexed atom-photon entanglement interface \- Optica Publishing Group, accessed March 22, 2025, [https://opg.optica.org/abstract.cfm?uri=oe-31-5-7200](https://opg.optica.org/abstract.cfm?uri=oe-31-5-7200)  
12. Extending the coherence time limit of a single-alkali-atom qubit by ..., accessed March 22, 2025, [https://opg.optica.org/optica/abstract.cfm?uri=optica-11-10-1391](https://opg.optica.org/optica/abstract.cfm?uri=optica-11-10-1391)  
13. Cavity-based quantum networks with single atoms and optical photons \- ResearchGate, accessed March 22, 2025, [https://www.researchgate.net/publication/269417316\_Cavity-based\_quantum\_networks\_with\_single\_atoms\_and\_optical\_photons](https://www.researchgate.net/publication/269417316_Cavity-based_quantum_networks_with_single_atoms_and_optical_photons)  
14. Cavity-Enhanced Optical Lattices for Scaling Neutral Atom Quantum Technologies to Higher Qubit Numbers \- Physical Review Link Manager, accessed March 22, 2025, [https://link.aps.org/doi/10.1103/PRXQuantum.3.030314](https://link.aps.org/doi/10.1103/PRXQuantum.3.030314)  
15. Noise-induced distributed entanglement in atom-cavity-fiber system, accessed March 22, 2025, [https://opg.optica.org/abstract.cfm?uri=oe-25-26-33359](https://opg.optica.org/abstract.cfm?uri=oe-25-26-33359)  
16. Fast delivery of heralded atom-photon quantum correlation over 12 km fiber through multiplexing enhancement \- PMC, accessed March 22, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11603145/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11603145/)  
17. Cavity QED with Atomic Ensembles \- Simon Lab, accessed March 22, 2025, [https://simonlab.stanford.edu/theses/SimonJon\_PhDThesis.pdf](https://simonlab.stanford.edu/theses/SimonJon_PhDThesis.pdf)  
18. Single-Atom Single-Photon Quantum Interface \- ResearchGate, accessed March 22, 2025, [https://www.researchgate.net/publication/6248381\_Single-Atom\_Single-Photon\_Quantum\_Interface](https://www.researchgate.net/publication/6248381_Single-Atom_Single-Photon_Quantum_Interface)  
19. Atom-mediated deterministic generation and stitching of photonic graph states \- arXiv, accessed March 22, 2025, [https://arxiv.org/html/2406.00860v3](https://arxiv.org/html/2406.00860v3)  
20. Research Details \- mpq.mpg.de, accessed March 22, 2025, [https://www.mpq.mpg.de/4996520/details](https://www.mpq.mpg.de/4996520/details)  
21. Scalable Networking of Neutral-Atom Qubits: Nanofiber-Based Approach for Multiprocessor Fault-Tolerant Quantum Computer, accessed March 22, 2025, [https://www.nano-qt.com/resources/20250204](https://www.nano-qt.com/resources/20250204)  
22. Error-Detected Quantum Operations with Neutral Atoms Mediated by an Optical Cavity, accessed March 22, 2025, [https://arxiv.org/html/2410.10787v1](https://arxiv.org/html/2410.10787v1)  
23. Scalable Networking of Neutral-Atom Qubits: Nanofiber-Based Approach for Multiprocessor Fault-Tolerant Quantum Computers \- Physical Review Link Manager, accessed March 22, 2025, [https://link.aps.org/doi/10.1103/PRXQuantum.6.010101](https://link.aps.org/doi/10.1103/PRXQuantum.6.010101)  
24. Error-Detected Quantum Operations with Neutral Atoms Mediated by an Optical Cavity, accessed March 22, 2025, [https://inspirehep.net/literature/2839898](https://inspirehep.net/literature/2839898)  
25. Generation and optimization of entanglement between atoms chirally coupled to spin cavities | Phys. Rev. Research \- Physical Review Link Manager, accessed March 22, 2025, [https://link.aps.org/doi/10.1103/PhysRevResearch.7.L012058](https://link.aps.org/doi/10.1103/PhysRevResearch.7.L012058)  
26. Tunable Ion-Photon Entanglement in an Optical Cavity \- PMC, accessed March 22, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC4337972/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4337972/)  
27. Deterministic spin-photon entanglement from a trapped ion in a fiber Fabry–Perot cavity \- Inspire HEP, accessed March 22, 2025, [https://inspirehep.net/files/13491e85769084b9190d8b0c0b7f0147](https://inspirehep.net/files/13491e85769084b9190d8b0c0b7f0147)  
28. Deterministic high-rate entanglement distillation with neutral atom arrays \- Inspire HEP, accessed March 22, 2025, [https://inspirehep.net/literature/2896282](https://inspirehep.net/literature/2896282)  
29. Deterministic spin-photon entanglement from a trapped ion in a fiber Fabry–Perot cavity, accessed March 22, 2025, [https://www.researchgate.net/publication/348697752\_Deterministic\_spin-photon\_entanglement\_from\_a\_trapped\_ion\_in\_a\_fiber\_Fabry-Perot\_cavity](https://www.researchgate.net/publication/348697752_Deterministic_spin-photon_entanglement_from_a_trapped_ion_in_a_fiber_Fabry-Perot_cavity)  
30. Shaping Quantum Noise through Cascaded Nonlinear Processes in a Dissipation-Engineered Multimode Cavity \- Physical Review Link Manager, accessed March 22, 2025, [https://link.aps.org/doi/10.1103/PRXQuantum.5.040345](https://link.aps.org/doi/10.1103/PRXQuantum.5.040345)  
31. Protecting the Quantum Coherence of Two Atoms Inside an Optical Cavity by Quantum Feedback Control Combined with Noise-Assisted Preparation \- MDPI, accessed March 22, 2025, [https://www.mdpi.com/2304-6732/11/5/400](https://www.mdpi.com/2304-6732/11/5/400)  
32. Cavity dark mode mediated by atom array without atomic scattering loss | Phys. Rev. Research \- Physical Review Link Manager, accessed March 22, 2025, [https://link.aps.org/doi/10.1103/PhysRevResearch.6.L042026](https://link.aps.org/doi/10.1103/PhysRevResearch.6.L042026)  
33. Cavity-enabled measurements and interactions in neutral atoms \- Elmore Family School of Electrical and Computer Engineering \- Purdue University, accessed March 22, 2025, [https://engineering.purdue.edu/ECE/Events/2025/cavity-enabled-measurements-and-interactions-in-neutral-atoms](https://engineering.purdue.edu/ECE/Events/2025/cavity-enabled-measurements-and-interactions-in-neutral-atoms)