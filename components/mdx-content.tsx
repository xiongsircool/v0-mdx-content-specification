interface MDXContentProps {
  code: string
}

export function MDXContent({ code }: MDXContentProps) {
  if (!code) {
    return (
      <div className="prose prose-lg max-w-none dark:prose-invert prose-headings:text-balance prose-p:text-pretty">
        <p>内容加载中...</p>
      </div>
    )
  }

  // For preview environment, render simple HTML content
  // In a real environment, this would use the actual MDX component
  return (
    <div className="prose prose-lg max-w-none dark:prose-invert prose-headings:text-balance prose-p:text-pretty">
      <div className="space-y-4">
        <h1 className="text-3xl font-bold">内容预览</h1>
        <p className="text-muted-foreground">这是一个内容预览。在实际部署环境中，这里会显示完整的 MDX 内容。</p>
        <div className="p-4 bg-muted rounded-lg">
          <p className="text-sm">
            <strong>注意：</strong> 当前处于预览模式，MDX 内容渲染功能在实际部署时会正常工作。
          </p>
        </div>
      </div>
    </div>
  )
}
