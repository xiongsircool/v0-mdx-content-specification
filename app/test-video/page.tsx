import { getFileBySlug } from "@/lib/server-markdown-loader"
import { MarkdownContent } from "@/components/markdown-content"
import { Card, CardContent } from "@/components/ui/card"

export default async function TestVideoPage() {
  // Load the video test post content
  const postFile = getFileBySlug('posts', '2025-10-06-video-link-test')
  const postContent = postFile?.htmlContent || '<p>内容未找到</p>'

  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto px-4 py-8 max-w-4xl">
        <h1 className="text-4xl font-bold mb-8">视频嵌入测试页面</h1>

        <Card>
          <CardContent className="p-8">
            {postContent ? (
              <MarkdownContent content={postContent} />
            ) : (
              <p>无法加载测试内容</p>
            )}
          </CardContent>
        </Card>
      </div>
    </div>
  )
}