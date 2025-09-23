"use client"

import type React from "react"

import { useMemo } from "react"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert"
import { Separator } from "@/components/ui/separator"
import { Info, AlertTriangle, CheckCircle, XCircle, Code, Copy, ExternalLink, FileText, Quote } from "lucide-react"
import { Button } from "@/components/ui/button"
import Image from "next/image"
import Link from "next/link"

interface MDXContentProps {
  code: string
  title?: string
}

// MDX 组件映射
const MDXComponents = {
  // 标题组件
  h1: ({ children, ...props }: any) => (
    <h1 className="text-3xl font-bold text-balance mb-6 mt-8 first:mt-0" {...props}>
      {children}
    </h1>
  ),
  h2: ({ children, ...props }: any) => (
    <h2 className="text-2xl font-semibold text-balance mb-4 mt-8 first:mt-0" {...props}>
      {children}
    </h2>
  ),
  h3: ({ children, ...props }: any) => (
    <h3 className="text-xl font-semibold text-balance mb-3 mt-6 first:mt-0" {...props}>
      {children}
    </h3>
  ),
  h4: ({ children, ...props }: any) => (
    <h4 className="text-lg font-semibold text-balance mb-2 mt-4 first:mt-0" {...props}>
      {children}
    </h4>
  ),

  // 段落和文本
  p: ({ children, ...props }: any) => (
    <p className="text-pretty leading-7 mb-4" {...props}>
      {children}
    </p>
  ),

  // 列表
  ul: ({ children, ...props }: any) => (
    <ul className="list-disc list-inside mb-4 space-y-2" {...props}>
      {children}
    </ul>
  ),
  ol: ({ children, ...props }: any) => (
    <ol className="list-decimal list-inside mb-4 space-y-2" {...props}>
      {children}
    </ol>
  ),
  li: ({ children, ...props }: any) => (
    <li className="leading-7" {...props}>
      {children}
    </li>
  ),

  // 链接
  a: ({ children, href, ...props }: any) => (
    <Link
      href={href}
      className="text-primary hover:text-primary/80 underline underline-offset-4 inline-flex items-center gap-1"
      {...props}
    >
      {children}
      {href?.startsWith("http") && <ExternalLink className="h-3 w-3" />}
    </Link>
  ),

  // 代码
  code: ({ children, className, ...props }: any) => {
    const isInline = !className
    if (isInline) {
      return (
        <code className="bg-muted px-2 py-1 rounded text-sm font-mono" {...props}>
          {children}
        </code>
      )
    }

    return (
      <Card className="mb-4">
        <CardHeader className="pb-3">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <Code className="h-4 w-4" />
              <span className="text-sm font-medium">代码示例</span>
            </div>
            <Button variant="ghost" size="sm" onClick={() => navigator.clipboard.writeText(children)}>
              <Copy className="h-4 w-4" />
            </Button>
          </div>
        </CardHeader>
        <CardContent className="pt-0">
          <pre className="bg-muted p-4 rounded-lg overflow-x-auto">
            <code className="text-sm font-mono" {...props}>
              {children}
            </code>
          </pre>
        </CardContent>
      </Card>
    )
  },

  // 引用
  blockquote: ({ children, ...props }: any) => (
    <Card className="mb-4 border-l-4 border-l-primary">
      <CardContent className="p-4">
        <div className="flex gap-3">
          <Quote className="h-5 w-5 text-primary flex-shrink-0 mt-1" />
          <div className="text-muted-foreground italic" {...props}>
            {children}
          </div>
        </div>
      </CardContent>
    </Card>
  ),

  // 分隔线
  hr: () => <Separator className="my-8" />,

  // 图片
  img: ({ src, alt, ...props }: any) => (
    <Card className="mb-4">
      <CardContent className="p-4">
        <div className="relative aspect-video overflow-hidden rounded-lg">
          <Image
            src={src || "/placeholder.svg"}
            alt={alt || "图片"}
            fill
            className="object-cover"
            unoptimized
            {...props}
          />
        </div>
        {alt && <p className="text-sm text-muted-foreground text-center mt-2">{alt}</p>}
      </CardContent>
    </Card>
  ),

  // 表格
  table: ({ children, ...props }: any) => (
    <Card className="mb-4 overflow-hidden">
      <CardContent className="p-0">
        <div className="overflow-x-auto">
          <table className="w-full" {...props}>
            {children}
          </table>
        </div>
      </CardContent>
    </Card>
  ),
  thead: ({ children, ...props }: any) => (
    <thead className="bg-muted" {...props}>
      {children}
    </thead>
  ),
  th: ({ children, ...props }: any) => (
    <th className="px-4 py-3 text-left font-medium" {...props}>
      {children}
    </th>
  ),
  td: ({ children, ...props }: any) => (
    <td className="px-4 py-3 border-t" {...props}>
      {children}
    </td>
  ),
}

// 自定义 MDX 组件
const CustomComponents = {
  // 提示框组件
  Callout: ({
    type = "info",
    title,
    children,
  }: { type?: "info" | "warning" | "success" | "error"; title?: string; children: React.ReactNode }) => {
    const config = {
      info: { icon: Info, className: "border-blue-200 bg-blue-50 text-blue-900" },
      warning: { icon: AlertTriangle, className: "border-yellow-200 bg-yellow-50 text-yellow-900" },
      success: { icon: CheckCircle, className: "border-green-200 bg-green-50 text-green-900" },
      error: { icon: XCircle, className: "border-red-200 bg-red-50 text-red-900" },
    }

    const { icon: Icon, className } = config[type]

    return (
      <Alert className={`mb-4 ${className}`}>
        <Icon className="h-4 w-4" />
        {title && <AlertTitle>{title}</AlertTitle>}
        <AlertDescription>{children}</AlertDescription>
      </Alert>
    )
  },

  // 步骤组件
  Steps: ({ children }: { children: React.ReactNode }) => (
    <div className="mb-6">
      <div className="space-y-4">{children}</div>
    </div>
  ),

  Step: ({ number, title, children }: { number: number; title: string; children: React.ReactNode }) => (
    <Card>
      <CardHeader className="pb-3">
        <CardTitle className="flex items-center gap-3">
          <div className="flex h-8 w-8 items-center justify-center rounded-full bg-primary text-primary-foreground text-sm font-bold">
            {number}
          </div>
          {title}
        </CardTitle>
      </CardHeader>
      <CardContent className="pt-0">{children}</CardContent>
    </Card>
  ),

  // 标签组件
  Tag: ({
    children,
    variant = "default",
  }: { children: React.ReactNode; variant?: "default" | "secondary" | "outline" }) => (
    <Badge variant={variant as any} className="mr-2 mb-2">
      {children}
    </Badge>
  ),

  // 文件组件
  File: ({ name, type, children }: { name: string; type?: string; children: React.ReactNode }) => (
    <Card className="mb-4">
      <CardHeader className="pb-3">
        <CardTitle className="flex items-center gap-2 text-base">
          <FileText className="h-4 w-4" />
          <code className="text-sm">{name}</code>
          {type && (
            <Badge variant="outline" className="text-xs">
              {type}
            </Badge>
          )}
        </CardTitle>
      </CardHeader>
      <CardContent className="pt-0">{children}</CardContent>
    </Card>
  ),
}

export function MDXContent({ code, title }: MDXContentProps) {
  const content = useMemo(() => {
    if (!code) {
      return (
        <div className="text-center py-8">
          <FileText className="h-12 w-12 mx-auto mb-4 text-muted-foreground opacity-50" />
          <p className="text-muted-foreground">内容加载中...</p>
        </div>
      )
    }

    // 在预览环境中，渲染示例内容
    return (
      <div className="space-y-6">
        <MDXComponents.h1>{title || "MDX 内容示例"}</MDXComponents.h1>

        <MDXComponents.p>
          这是一个 MDX 内容渲染示例。在实际部署环境中，这里会显示完整的 MDX 内容，支持丰富的组件和交互功能。
        </MDXComponents.p>

        <CustomComponents.Callout type="info" title="预览模式">
          当前处于预览模式，展示 MDX 组件的渲染效果。实际部署时会正确解析和渲染 MDX 文件内容。
        </CustomComponents.Callout>

        <MDXComponents.h2>功能特性</MDXComponents.h2>

        <MDXComponents.ul>
          <MDXComponents.li>支持标准 Markdown 语法</MDXComponents.li>
          <MDXComponents.li>自定义 React 组件集成</MDXComponents.li>
          <MDXComponents.li>代码高亮和复制功能</MDXComponents.li>
          <MDXComponents.li>响应式图片和媒体</MDXComponents.li>
        </MDXComponents.ul>

        <MDXComponents.h3>代码示例</MDXComponents.h3>

        <MDXComponents.code className="language-javascript">
          {`function HelloWorld() {
  return (
    <div>
      <h1>Hello, MDX!</h1>
      <p>这是一个 React 组件示例</p>
    </div>
  )
}`}
        </MDXComponents.code>

        <CustomComponents.Steps>
          <CustomComponents.Step number={1} title="安装依赖">
            <MDXComponents.p>首先安装必要的依赖包</MDXComponents.p>
            <MDXComponents.code>npm install @next/mdx</MDXComponents.code>
          </CustomComponents.Step>

          <CustomComponents.Step number={2} title="配置 MDX">
            <MDXComponents.p>在 next.config.js 中配置 MDX 支持</MDXComponents.p>
          </CustomComponents.Step>

          <CustomComponents.Step number={3} title="创建内容">
            <MDXComponents.p>开始创建你的 MDX 内容文件</MDXComponents.p>
          </CustomComponents.Step>
        </CustomComponents.Steps>

        <MDXComponents.blockquote>
          MDX 让你能够在 Markdown 中无缝使用 JSX，为技术文档和博客内容提供了强大的表达能力。
        </MDXComponents.blockquote>

        <CustomComponents.Callout type="success" title="部署就绪">
          所有 MDX 组件已配置完成，可以在生产环境中正常使用。
        </CustomComponents.Callout>
      </div>
    )
  }, [code, title])

  return (
    <div className="prose prose-lg max-w-none dark:prose-invert prose-headings:text-balance prose-p:text-pretty">
      {content}
    </div>
  )
}
