# Meruki_CD_ShippingCostOptimizer
为挖煤姬购买CD专辑设计的合单方案优化器，最小化运费成本

<img width="2560" height="1520" alt="image" src="https://github.com/user-attachments/assets/9a277127-b7b8-41aa-a404-b7cd18f603ef" />

---

程序专为海淘CD设计，不建议将其他物品与CD混发，可能增加被税风险。

程序支持5kg以下的合单包裹，支持的运输方式如下：

<img width="1742" height="966" alt="image" src="https://github.com/user-attachments/assets/f18454a1-71ff-4797-8d88-bd44be1ecf91" />

程序不支持使用竹蜻蜓PLUS发送体积大于9000的包裹.

---

在运行优化器前，你需要正确填写`ShoppingList.csv`，一行对应一个待合单订单，示例如下：

<img width="1701" height="779" alt="image" src="https://github.com/user-attachments/assets/0a98d9f7-0d0f-4f82-a17c-a8500ca2fecc" />

根据海关规定，单碟（盘）发行的音像制品，每人每次 20 盘以下免税，但是没有解释类似“初回盘”的 CD+DVD 专辑是按照一盘还是两盘计算。

保险起见，建议将此类 CD 作为单碟（盘）发行的音像制品，按照 2 盘计算。即“专辑数”+1，“盘数”+2。

请注意，成套发行的音像制品，每人每次 3 套以下免税。对于含有多张 CD 碟片的专辑，有被认为是成套音像制品的风险，谨慎使用此程序。

---

`run_pso_sa_optimization`为`false`时，仅使用贪心算法求解；为`true`时，程序会在贪心结果的基础上使用 PSO-SA 进一步优化。

本程序支持积分抵扣计算，请在`points_balance`填入你的积分余额。

---

**⚠️⚠️注意⚠️⚠️**

程序没有对一些物流方式的细节限制做出判断

(譬如：竹蜻蜓PLUS要求包裹长度不超过60cm；竹蜻蜓Max要求包裹金额不得高于2000RMB等。)

请务必仔细阅读官方使用指南: https://www.meruki.cn/userGuide/novice-strategy?code=novice-strategy ，确保合单方案合规。

由于平台规则、海关政策可能发生改变，此程序的计算结果仅供参考。
