# 🧬 SBC - 上海生物信息中心学生组织官网

<div align="center">

![SBC Logo](public/sbc-logo.svg)

**上海生物信息中心学生组织官方网站**

一个现代化的生物信息学学习社区平台，为学生提供技术分享、学习资源和活动信息。

[![Next.js](https://img.shields.io/badge/Next.js-15-black?style=for-the-badge&logo=next.js)](https://nextjs.org/)
[![TypeScript](https://img.shields.io/badge/TypeScript-5-blue?style=for-the-badge&logo=typescript)](https://www.typescriptlang.org/)
[![Tailwind CSS](https://img.shields.io/badge/Tailwind%20CSS-v4-cyan?style=for-the-badge&logo=tailwindcss)](https://tailwindcss.com/)
[![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)](LICENSE)

</div>

---

## 🌟 项目简介

SBC（上海生物信息中心）学生组织官网是一个基于 Next.js 15 构建的现代化内容管理平台，专门为生物信息学学习社区设计。网站采用自定义 Markdown 内容管理系统，支持技术分享、活动公告、学习资源和会议信息管理。

### ✨ 主要特性

- 🎨 **现代化设计**：基于 Tailwind CSS v4 和 shadcn/ui 组件库
- 📱 **响应式布局**：完美适配桌面端和移动端
- 🧬 **DNA主题设计**：融合生物信息学特色的视觉元素
- 🔍 **全局搜索**：支持跨内容类型的实时搜索功能
- 📊 **智能统计**：动态展示社区活动和内容统计
- 🗓️ **会议管理**：自动状态计算和时间排序
- 🎯 **内容分类**：公告、技术分享、学习资料、会议信息
- 🌓 **暗色模式**：支持明暗主题切换
- ⚡ **高性能**：SSG + SSR 混合渲染模式

---

## 🖼️ 网站效果预览

### 网站实际演示
![网站效果预览](images/网站实际效果预览.gif)

### 功能特性展示
- **响应式设计**: 完美适配桌面端和移动端
- **现代UI**: 基于Tailwind CSS v4的美观界面
- **DNA主题设计**: 生物信息学特色的视觉元素
- **交互动画**: Framer Motion驱动的流畅体验
- **全局搜索**: 支持跨内容类型的实时搜索
- **智能统计**: 动态展示社区活动和内容统计
- **会议管理**: 自动状态计算和时间排序

### 设计亮点
- **渐变配色**: 蓝色→紫色→粉色的现代渐变
- **组件库**: shadcn/ui + Radix UI的无障碍设计
- **代码高亮**: highlight.js驱动的语法高亮
- **Logo设计**: DNA双螺旋结构，体现生物信息特色

---

## 🎯 网站功能

### 1. 首页展示
- **实时统计数据**：活跃成员、技术分享、学习资源、社区活动
- **最新公告**：优先级排序的过期公告自动过滤
- **技术分享**：最新技术文章和项目成果展示
- **文献展示**：外部生物信息学文献资源
- **合作伙伴**：赞助商和合作机构展示

### 2. 内容管理
#### 📢 活动公告
- 招新信息、比赛通知、讲座预告、活动安排
- 自动过期机制和优先级排序
- 支持附件下载和外部链接

#### 📝 技术分享
- 多作者支持的技术博客文章
- 代码高亮和复制功能
- 标签分类和筛选功能
- 封面图片和摘要展示

#### 📚 学习资料
- 文章、视频、项目、数据集、幻灯片分类
- 入门、进阶、实战难度等级
- 学习路径和先修要求
- 外部资源链接管理

#### 🗓️ 会议信息
- 自动状态计算（即将开始/进行中/已结束）
- 时间优先级排序
- 在线会议链接和容量管理
- 材料下载和报名统计

### 3. 交互功能
- **全局搜索**：支持标题、摘要、标签、作者搜索
- **智能筛选**：按状态、类型、标签等多维度筛选
- **响应式导航**：移动端优化的折叠菜单
- **动画效果**：流畅的页面过渡和交互动画

---

## 🛠️ 技术栈

### 前端框架
- **Next.js 15** - React 全栈框架，支持 App Router
- **React 19** - 用户界面构建库
- **TypeScript** - 类型安全的 JavaScript 超集

### 样式和组件
- **Tailwind CSS v4** - 实用优先的 CSS 框架
- **shadcn/ui** - 现代化的 React 组件库
- **Radix UI** - 无障碍访问的组件基础
- **Framer Motion** - 动画和手势库

### 内容管理
- **gray-matter** - Markdown 前置元数据解析
- **marked** - Markdown 转 HTML 解析器
- **highlight.js** - 代码语法高亮
- **自定义加载系统** - 服务端 Markdown 文件处理

### 开发工具
- **pnpm** - 快速、节省磁盘空间的包管理器
- **ESLint** - JavaScript 代码质量检查
- **PostCSS** - CSS 后处理器
- **Autoprefixer** - CSS 前缀自动添加

### 部署和分析
- **Vercel** - 现代化应用部署平台
- **Vercel Analytics** - 网站访问分析
- **GitHub Actions** - CI/CD 自动化部署

---

## 📦 项目结构

```
SBCBIOchina/
├── app/                     # Next.js App Router 页面
│   ├── about/               # 关于我们页面
│   ├── admin/               # 管理后台
│   ├── announcements/       # 公告详情页
│   ├── api/                 # API 路由
│   ├── learn/               # 学习资源页面
│   ├── meetings/            # 会议页面
│   ├── posts/               # 技术分享页面
│   ├── search/              # 搜索页面
│   ├── [slug]/              # 动态路由
│   ├── layout.tsx           # 根布局
│   └── page.tsx             # 首页
├── components/              # React 组件
│   ├── animations.tsx       # 动画组件
│   ├── navigation.tsx       # 导航组件
│   ├── content-card.tsx     # 内容卡片
│   ├── features-grid.tsx    # 功能展示
│   └── ...                  # 其他组件
├── lib/                     # 工具库
│   ├── server-markdown-loader.ts  # Markdown 加载器
│   └── post-manager.ts      # 文章管理器
├── content/                 # Markdown 内容
│   ├── announcements/       # 公告内容
│   ├── posts/              # 技术分享
│   ├── learn/              # 学习资源
│   └── meetings/           # 会议信息
├── public/                 # 静态资源
│   ├── logo-nav.svg        # 导航 Logo
│   ├── sbc-logo.svg        # 主 Logo
│   └── ...                 # 其他静态文件
├── middleware.ts           # 中间件
├── tailwind.config.ts      # Tailwind 配置
├── next.config.mjs         # Next.js 配置
└── package.json            # 项目依赖
```

---

## 🚀 快速开始

### 环境要求

- **Node.js**: >= 18.0.0
- **pnpm**: >= 8.0.0
- **Git**: 最新版本

### 安装步骤

1. **克隆项目**
```bash
git clone https://github.com/your-username/SBCBIOchina.git
cd SBCBIOchina
```

2. **安装依赖**
```bash
pnpm install
```

3. **启动开发服务器**
```bash
pnpm dev
```

4. **打开浏览器**
访问 [http://localhost:3000](http://localhost:3000) 查看网站。

### 开发脚本

```bash
# 开发模式
pnpm dev

# 构建生产版本
pnpm build

# 启动生产服务器
pnpm start

# 代码检查
pnpm lint

# 类型检查
npx tsc --noEmit
```

---

## 📝 内容管理指南

### 内容结构

所有内容都存储在 `content/` 目录下的 Markdown 文件中，使用 YAML 前置元数据：

#### 公告示例
```markdown
---
title: "2025年秋季招新公告"
publishedAt: "2025-09-28"
priority: 1
tags: ["招新"]
expiresAt: "2025-10-31"
summary: "上海生物信息中心2025年秋季招新开始，欢迎对生物信息学感兴趣的同学加入我们。"
---

# 招新公告内容...

## 报名方式
- **报名时间**: 2025年9月28日 - 10月31日
- **报名方式**: 在线报名表单
```

#### 技术分享示例
```markdown
---
title: "单细胞RNA测序数据分析指南"
publishedAt: "2025-10-05"
updatedAt: "2025-10-05"
tags: ["单细胞", "RNA-seq", "数据分析"]
authors: ["张三", "李四"]
excerpt: "详细介绍单细胞RNA测序数据的分析流程和常用工具"
coverImage: "/images/single-cell-cover.jpg"
---

# 单细胞RNA测序数据分析指南

## 分析流程
1. 数据质控
2. 细胞聚类
3. 差异表达分析
```

### 自动化功能

- **会议状态**: 基于时间自动计算（即将开始/进行中/已结束）
- **公告过期**: 超过 `expiresAt` 日期自动标记为过期
- **URL生成**: 根据文件路径自动生成页面URL
- **代码高亮**: 自动语法高亮和复制按钮

---

## 🎨 设计系统

### 颜色主题

```css
/* 主色调 */
--primary: #3B82F6;    /* 蓝色 */
--secondary: #8B5CF6;  /* 紫色 */
--accent: #EC4899;     /* 粉色 */

/* 渐变组合 */
primary-gradient: linear-gradient(135deg, #3B82F6, #8B5CF6);
secondary-gradient: linear-gradient(135deg, #8B5CF6, #EC4899);
```

### 组件库

项目使用以下组件库：
- **shadcn/ui**: 基础 UI 组件
- **Radix UI**: 无头组件基础
- **Lucide React**: 图标库
- **Framer Motion**: 动画效果

### Logo 设计

![Logo 设计说明](public/sbc-logo.svg)

- **DNA双螺旋**: 生物信息学的核心象征
- **渐变配色**: 体现科技感和现代性
- **连接结构**: 代表数据连接和协作
- **圆形背景**: 象征完整性和社区包容性

---

## 🚀 部署指南

### Vercel 部署（推荐）

1. **推送代码到 GitHub**
```bash
git add .
git commit -m "feat: add new features"
git push origin main
```

2. **连接 Vercel**
   - 登录 [Vercel](https://vercel.com)
   - 导入 GitHub 仓库
   - 自动部署

3. **环境变量设置**
   - 在 Vercel 控制台设置环境变量
   - 配置分析工具等第三方服务

### 其他平台部署

```bash
# 构建项目
pnpm build

# 导出静态文件（如需）
npx next export

# 启动服务
pnpm start
```

---

## 🤝 贡献指南

### 开发流程

1. **Fork 项目**
2. **创建功能分支**
   ```bash
   git checkout -b feature/new-feature
   ```

3. **提交更改**
   ```bash
   git commit -m "feat: add new feature"
   ```

4. **推送分支**
   ```bash
   git push origin feature/new-feature
   ```

5. **创建 Pull Request**

### 代码规范

- 使用 TypeScript 进行类型定义
- 遵循 ESLint 规则
- 组件使用函数式写法
- 样式使用 Tailwind CSS
- 提交信息符合 Conventional Commits

### 内容贡献

- 在 `content/` 目录下创建 Markdown 文件
- 遵循现有文件命名规范
- 完善前置元数据
- 确保内容质量

---

## 📈 性能优化

### 已实现的优化

- **图片优化**: Next.js 自动图片优化
- **代码分割**: 路由级别的代码分割
- **缓存策略**: 静态资源缓存
- **SEO 优化**: 元标签和结构化数据
- **懒加载**: 图片和组件懒加载

### 性能指标

- **Lighthouse Score**: >90
- **Core Web Vitals**: 全部通过
- **加载速度**: <2s
- **首次绘制**: <1s

---

## 🔧 自定义配置

### 修改网站信息

编辑 `app/layout.tsx` 中的元数据：

```typescript
export const metadata: Metadata = {
  title: "SBC - 上海生物信息中心学生组织",
  description: "上海生物信息中心学生组织官方网站...",
  // 其他元数据...
}
```

### 添加新内容类型

1. 在 `lib/server-markdown-loader.ts` 中添加加载函数
2. 创建对应的页面路由
3. 更新搜索和筛选功能

### 主题定制

修改 `tailwind.config.ts` 配置文件，或更新 CSS 变量。

---

## 🐛 故障排除

### 常见问题

**Q: 依赖安装失败**
```bash
# 清除缓存重新安装
rm -rf node_modules package-lock.json
pnpm install
```

**Q: 开发服务器启动失败**
```bash
# 检查端口占用
lsof -i :3000

# 使用其他端口
pnpm dev -p 3001
```

**Q: 构建失败**
```bash
# 检查 TypeScript 错误
npx tsc --noEmit

# 检查 ESLint 错误
pnpm lint
```

### 调试技巧

- 使用浏览器开发者工具检查网络请求
- 查看 Next.js 构建日志
- 检查中间件配置
- 验证 Markdown 文件格式

---

## 📄 许可证

本项目采用 MIT 许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

---

## 🙏 致谢

感谢以下开源项目和社区的支持：

- [Next.js](https://nextjs.org/) - React 框架
- [Tailwind CSS](https://tailwindcss.com/) - CSS 框架
- [shadcn/ui](https://ui.shadcn.com/) - 组件库
- [Vercel](https://vercel.com/) - 部署平台

特别感谢所有为 SBC 学生组织贡献内容的成员！

---

## 📞 联系我们

- **官方网站**: [https://sbcbiochina.org](https://sbcbiochina.org)
- **邮箱**: [contact@sbcbiochina.org](mailto:contact@sbcbiochina.org)
- **GitHub**: [SBCBIOchina](https://github.com/your-username/SBCBIOchina)

<div align="center">

![SBC Logo](public/logo-nav.svg)

**Made with ❤️ by SBC Student Organization**

</div>
